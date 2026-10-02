"""
Run one generative model on one ChEMBL seed split.

Takes results/ChEMBL_splits/split_<NNN>.csv (from 01_prepare_chembl_splits.py)
and runs the model's own run.sh on it, inside the CPU env built by
03_build_cpu_envs.py, writing
--path-to-output/<model-id>/csv/split_<NNN>.csv (split_<NNN>_repeated.csv with
--repeated), and, under logs/ in the same model folder, a .log with the model's
stdout/stderr and a .json with the run's provenance (wall time, host, CPUs,
Slurm ids, model commit, exit code). Slurm's own .out/.err files go to out/
(see 05_submit_slurm.py).

The env is used directly (its bin/ first on PATH), not through `conda run`, so
the script also works on machines without conda, e.g. cluster nodes. The user
site is disabled (PYTHONNOUSERSITE=1) so the env cannot lean on the caller's
~/.local, and if --hf-home exists it is used as an offline HuggingFace cache
(HF_HOME + HF_HUB_OFFLINE=1). A finished output is never overwritten unless
--force, so a job array can be re-submitted safely. The exit code is the
model's, or 1 if its output does not have one row per input.

Standard library only; fully self-contained (no chemsampler imports). Runs on
Python 3.8+, so it can be launched with any model env's interpreter.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import socket
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AUDIT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DEFAULT_MODELS_DIR = os.path.join(AUDIT_DIR, "models")
DEFAULT_ENVS_DIR = os.path.join(AUDIT_DIR, "envs_cpu")
DEFAULT_SPLITS_DIR = os.path.join(AUDIT_DIR, "results", "ChEMBL_splits")
DEFAULT_OUTPUT_DIR = os.path.join(AUDIT_DIR, "results")
DEFAULT_HF_HOME = os.path.join(AUDIT_DIR, "cache", "huggingface")


def available_cpus() -> int:
    """CPUs this process may use: the Slurm allocation if any, else the
    affinity mask (which respects cpusets), else every core of the machine."""
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm:
        return int(slurm)
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def count_rows(path: str) -> int:
    """Number of data rows (header excluded) in a CSV file."""
    with open(path, newline="") as f:
        return max(sum(1 for _ in csv.reader(f)) - 1, 0)


def model_commit(model_dir: str) -> str | None:
    """Short git commit of the model clone, or None if it cannot be read."""
    result = subprocess.run(
        ["git", "-C", model_dir, "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() or None if result.returncode == 0 else None


def build_env(env_dir: str, cpus: int, hf_home: str) -> dict[str, str]:
    """Environment for the model process: its env first on PATH, user site off,
    BLAS/OpenMP threads pinned to the allotted CPUs, offline HF cache if any."""
    env = {**os.environ}
    env["PATH"] = os.path.join(env_dir, "bin") + os.pathsep + env["PATH"]
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    env["OMP_NUM_THREADS"] = env["MKL_NUM_THREADS"] = str(cpus)
    if os.path.isdir(hf_home):
        env["HF_HOME"] = hf_home
        env["HF_HUB_OFFLINE"] = "1"
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Model identifier, e.g. eos2401")
    parser.add_argument(
        "--split-id", type=int, required=True, help="Split number (1-100)"
    )
    parser.add_argument(
        "--repeated",
        action="store_true",
        help="Save as split_<NNN>_repeated.csv (a second run of the same split)",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=None,
        help="Threads for the model (default: Slurm allocation, else all "
        "CPUs this process may use)",
    )
    parser.add_argument(
        "--force", action="store_true", help="Re-run even if the output exists"
    )
    parser.add_argument("--path-to-models", default=DEFAULT_MODELS_DIR)
    parser.add_argument("--path-to-envs", default=DEFAULT_ENVS_DIR)
    parser.add_argument("--path-to-splits", default=DEFAULT_SPLITS_DIR)
    parser.add_argument("--path-to-output", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--hf-home", default=DEFAULT_HF_HOME)
    args = parser.parse_args()

    model_dir = os.path.join(args.path_to_models, args.model)
    framework_dir = os.path.join(model_dir, "model", "framework")
    env_dir = os.path.join(args.path_to_envs, f"{args.model}-cpu")
    split_name = f"split_{args.split_id:03d}"
    split_path = os.path.join(args.path_to_splits, f"{split_name}.csv")
    stem = split_name + ("_repeated" if args.repeated else "")
    model_out = os.path.join(args.path_to_output, args.model)
    csv_dir = os.path.join(model_out, "csv")
    logs_dir = os.path.join(model_out, "logs")
    output = os.path.join(csv_dir, f"{stem}.csv")
    partial = os.path.join(csv_dir, f"{stem}.part.csv")

    for label, path in (
        ("run.sh", os.path.join(framework_dir, "run.sh")),
        ("env", os.path.join(env_dir, "bin", "python")),
        ("split", split_path),
    ):
        if not os.path.exists(path):
            raise SystemExit(f"{label} not found: {path}")

    if not args.force and os.path.exists(output) and os.path.getsize(output) > 0:
        print(f"{args.model} {stem}: output exists, skipping (use --force)")
        return

    os.makedirs(csv_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    cpus = args.cpus or available_cpus()
    n_inputs = count_rows(split_path)
    print(
        f"{args.model} {stem}: {n_inputs} inputs, {cpus} CPUs, host {socket.gethostname()}"
    )

    started = datetime.datetime.now(datetime.timezone.utc)
    t0 = time.monotonic()
    with open(os.path.join(logs_dir, f"{stem}.log"), "w") as log:
        result = subprocess.run(
            [
                "bash",
                os.path.join(framework_dir, "run.sh"),
                framework_dir,
                split_path,
                partial,
            ],
            env=build_env(env_dir, cpus, args.hf_home),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    wall = time.monotonic() - t0

    exit_code = result.returncode
    n_outputs = count_rows(partial) if os.path.exists(partial) else None
    if exit_code == 0 and n_outputs != n_inputs:
        print(f"  output has {n_outputs} rows for {n_inputs} inputs")
        exit_code = 1
    if exit_code == 0:
        os.replace(partial, output)
    elif os.path.exists(partial):
        os.remove(partial)

    provenance = {
        "model": args.model,
        "model_commit": model_commit(model_dir),
        "split": split_name,
        "repeated": args.repeated,
        "n_inputs": n_inputs,
        "n_output_rows": n_outputs,
        "exit_code": exit_code,
        "wall_seconds": round(wall, 1),
        "started_utc": started.isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "cpus": cpus,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "env": env_dir,
    }
    with open(os.path.join(logs_dir, f"{stem}.json"), "w") as f:
        json.dump(provenance, f, indent=2)
        f.write("\n")

    if exit_code == 0:
        print(f"  ok in {wall:.0f}s -> {output}")
    else:
        print(f"  FAILED (exit {exit_code}) in {wall:.0f}s, see logs/{stem}.log")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
