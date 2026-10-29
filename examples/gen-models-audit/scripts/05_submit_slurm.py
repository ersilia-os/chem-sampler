"""
Submit a Slurm job array that runs one model over every ChEMBL seed split.

Each array task calls 04_run_split.py with the model's own env interpreter
(no conda, no system Python needed on the node): task N (1..n_splits) runs
split N, and one extra last task re-runs --repeat-split as
split_<NNN>_repeated.csv. Everything lands in --path-to-output/<model-id>/:
csv/ (outputs), logs/ (model logs, provenance .json, and the batch scripts
submitted) and out/ (Slurm's own .out and .err files); see 04_run_split.py.

Only the lab's own nodes may be used: other nodes can be billed to the PI.
They are not written in this file: set CHEMSAMPLER_LAB_CPU_NODES and
CHEMSAMPLER_LAB_GPU_NODES (comma-separated node names, e.g. in your shell
profile). The script refuses any --nodelist entry outside those lists, and
refuses to run when they are not set, unless --allow-any-node is given. The
default --nodelist is the CPU list, on the preemptible spot_cpu partition, so a
job yields to other users instead of blocking them (tasks are requeued). The
nodelist + --nodes=1 form is the one the lab's other array scripts use.

sbatch must be reachable: run this on the cluster, or from any machine that
shares the filesystem with it, in which case the submission goes through
`ssh` to CHEMSAMPLER_SSH_HOST (or --ssh-host). --dry-run prints the batch
script and changes nothing.

Standard library only; fully self-contained (no chemsampler imports).
"""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import re
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AUDIT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
RUN_SPLIT = os.path.join(SCRIPT_DIR, "04_run_split.py")
DEFAULT_ENVS_DIR = os.path.join(AUDIT_DIR, "envs_cpu")
DEFAULT_SPLITS_DIR = os.path.join(AUDIT_DIR, "results", "ChEMBL_splits")
DEFAULT_OUTPUT_DIR = os.path.join(AUDIT_DIR, "results")

# The lab's own nodes (free for us); everything else is billed to the PI.


def env_list(name: str) -> list[str]:
    """Comma-separated environment variable as a list (empty if unset)."""
    return [x.strip() for x in os.environ.get(name, "").split(",") if x.strip()]


LAB_CPU_NODES = env_list("CHEMSAMPLER_LAB_CPU_NODES")
LAB_GPU_NODES = env_list("CHEMSAMPLER_LAB_GPU_NODES")
LAB_NODES = LAB_CPU_NODES + LAB_GPU_NODES
DEFAULT_SSH_HOST = os.environ.get("CHEMSAMPLER_SSH_HOST")
ARRAY_RE = re.compile(r"^[0-9][0-9,\-]*(%[0-9]+)?$")

SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={model}-bench
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition={partition}
#SBATCH --nodelist={nodelist}
#SBATCH --array={array}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --requeue
#SBATCH --output={slurm_dir}/%x_%A_%a.out
#SBATCH --error={slurm_dir}/%x_%A_%a.err
set -euo pipefail
cd "{audit_dir}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1

N_SPLITS={n_splits}
REPEAT_SPLIT={repeat_split}

task=$SLURM_ARRAY_TASK_ID
if [ "$task" -le "$N_SPLITS" ]; then
  which="--split-id $task"
else
  which="--split-id $REPEAT_SPLIT --repeated"
fi

"{env_python}" "{run_split}" --model {model} $which --cpus {cpus} \\
  --path-to-envs "{envs_dir}" --path-to-splits "{splits_dir}" \\
  --path-to-output "{output_dir}"
"""


def check_lab_only(nodelist: list[str]) -> None:
    """Exit with an error unless the request stays on lab-owned nodes."""
    if not LAB_NODES:
        raise SystemExit(
            "the lab's nodes are not configured: set CHEMSAMPLER_LAB_CPU_NODES "
            "and CHEMSAMPLER_LAB_GPU_NODES (comma-separated), or pass --allow-any-node"
        )
    if not nodelist:
        raise SystemExit(f"--nodelist is required (lab nodes: {','.join(LAB_NODES)})")
    outside = [n for n in nodelist if n not in LAB_NODES]
    if outside:
        raise SystemExit(
            f"not lab nodes: {', '.join(outside)} (lab nodes: {','.join(LAB_NODES)}); "
            "pass --allow-any-node to override"
        )


def submit(script_path: str, ssh_host: str | None) -> str:
    """sbatch the script (over ssh when sbatch is not local); return the job id."""
    if shutil.which("sbatch"):
        cmd = ["sbatch", "--parsable", script_path]
    else:
        if not ssh_host:
            raise SystemExit(
                "sbatch is not available here: set CHEMSAMPLER_SSH_HOST or pass --ssh-host"
            )
        cmd = [
            "ssh",
            "-o",
            "BatchMode=yes",
            ssh_host,
            "sbatch",
            "--parsable",
            script_path,
        ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise SystemExit(f"submission failed: {result.stderr.strip()}")
    return result.stdout.strip().split(";")[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Model identifier, e.g. eos2401")
    parser.add_argument("--partition", default="spot_cpu")
    parser.add_argument(
        "--nodelist",
        default=",".join(LAB_CPU_NODES),
        help="Comma-separated nodes (default: $CHEMSAMPLER_LAB_CPU_NODES)",
    )
    parser.add_argument("--cpus", type=int, default=4, help="CPUs per task")
    parser.add_argument("--mem", default="8G", help="Memory per task")
    parser.add_argument("--time", default="00:30:00", help="Time limit per task")
    parser.add_argument(
        "--repeat-split",
        type=int,
        default=1,
        help="Split that gets one extra run saved as _repeated (0 = none)",
    )
    parser.add_argument(
        "--array",
        default=None,
        help="Slurm array spec, e.g. 1 for a smoke test (default: every task, "
        "throttled by --max-parallel)",
    )
    parser.add_argument(
        "--max-parallel", type=int, default=50, help="Tasks running at once"
    )
    parser.add_argument(
        "--ssh-host",
        default=DEFAULT_SSH_HOST,
        help="Host to submit through when sbatch is not local (default: $CHEMSAMPLER_SSH_HOST)",
    )
    parser.add_argument("--path-to-envs", default=DEFAULT_ENVS_DIR)
    parser.add_argument("--path-to-splits", default=DEFAULT_SPLITS_DIR)
    parser.add_argument("--path-to-output", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--allow-any-node",
        action="store_true",
        help="Skip the lab-nodes-only check (nodes outside the lab cost money)",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print the batch script, submit nothing"
    )
    args = parser.parse_args()

    nodelist = [n for n in args.nodelist.split(",") if n]
    if not args.allow_any_node:
        check_lab_only(nodelist)

    env_python = os.path.join(args.path_to_envs, f"{args.model}-cpu", "bin", "python")
    if not os.path.exists(env_python):
        raise SystemExit(f"env not found: {env_python}")
    n_splits = len(glob.glob(os.path.join(args.path_to_splits, "split_*.csv")))
    if n_splits == 0:
        raise SystemExit(f"no split_*.csv in {args.path_to_splits}")
    if not 0 <= args.repeat_split <= n_splits:
        raise SystemExit(f"--repeat-split must be 0..{n_splits}")

    n_tasks = n_splits + (1 if args.repeat_split else 0)
    array = args.array or f"1-{n_tasks}%{args.max_parallel}"
    if not ARRAY_RE.match(array):
        raise SystemExit(f"invalid --array {array!r}")

    model_out = os.path.join(args.path_to_output, args.model)
    logs_dir = os.path.join(model_out, "logs")
    slurm_dir = os.path.join(model_out, "out")
    script = SBATCH_TEMPLATE.format(
        model=args.model,
        partition=args.partition,
        nodelist=",".join(nodelist),
        array=array,
        cpus=args.cpus,
        mem=args.mem,
        time=args.time,
        slurm_dir=slurm_dir,
        audit_dir=AUDIT_DIR,
        n_splits=n_splits,
        repeat_split=args.repeat_split or n_splits + 1,
        env_python=env_python,
        run_split=RUN_SPLIT,
        envs_dir=args.path_to_envs,
        splits_dir=args.path_to_splits,
        output_dir=args.path_to_output,
    )

    if args.dry_run:
        print(script)
        print(
            f"# dry run: {n_splits} splits + {n_tasks - n_splits} repeat, nothing written or submitted"
        )
        return

    os.makedirs(slurm_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d-%H%M%S")
    script_path = os.path.join(logs_dir, f"submit_{stamp}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)

    job_id = submit(script_path, args.ssh_host)
    print(f"submitted job {job_id} (array {array}) from {script_path}")
    print(f"  watch:  ssh {args.ssh_host} squeue -j {job_id}")
    print(
        f"  after:  ssh {args.ssh_host} sacct -j {job_id} --format=JobID,State,Elapsed,NodeList,ExitCode"
    )
    print(
        f"  slurm .out/.err: {slurm_dir}/   model logs: {logs_dir}/   csv: {model_out}/csv/"
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
