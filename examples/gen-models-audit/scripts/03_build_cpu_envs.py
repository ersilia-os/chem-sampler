"""
Build a CPU-only conda environment for each cloned generative model.

For every model directory under --path-to-models (as populated by
02_fetch_generative_models.py), builds a conda env at
--path-to-envs/<model-id>-cpu from that model's own declared dependencies.
Two dependency-spec formats exist across the Ersilia Model Hub:

- install.yml (newer models): an explicit "python" version plus a
  "commands" list of ["pip"|"conda", name, version?, *extra_args] entries.
- Dockerfile (older models): Python version comes from the base image tag
  (e.g. bentoml/model-server:0.11.0-py310 -> 3.10); every RUN instruction
  is executed literally inside the new env, since these lines are not all
  simple package installs (e.g. eos4qda pipes a shell script through
  `curl | bash`).

Every build step runs with PYTHONNOUSERSITE=1. Without it pip sees the
packages in the builder's ~/.local as "already satisfied" and skips them, so
the env only works on machines that share that home (it breaks on a cluster
node). After the steps, `pip check` (user site off) warns about requirements
that are still not installed.

Requires `conda`, `git` (only for cloning, done by script 02) and PyYAML
(the `audit` extra: `pip install -e ".[audit]"`) on PATH/importable. Fully
self-contained: no imports from chemsampler.
"""

import argparse
import os
import re
import shutil
import subprocess

import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODELS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "models"))
DEFAULT_ENVS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "envs_cpu"))

DOCKERFILE_PY_RE = re.compile(r"-py(\d)(\d+)")
CONDA_INSTALL_RE = re.compile(r"^\s*conda\s+install\b")
CONDA_YES_RE = re.compile(r"(^|\s)(-y|--yes)(\s|$)")

# Keep the builder's ~/.local out of every step, see the module docstring.
ISOLATED_ENV = {**os.environ, "PYTHONNOUSERSITE": "1"}


def find_models(models_dir: str, only: list[str] | None) -> list[str]:
    """Return model identifiers to process: `only` if given, else every
    subdirectory of `models_dir`."""
    if only:
        return only
    return sorted(
        name
        for name in os.listdir(models_dir)
        if os.path.isdir(os.path.join(models_dir, name))
    )


def load_install_yml(model_dir: str) -> tuple[str, list[list[str]]]:
    """Parse install.yml into (python_version, [argv, ...]) for `conda run`."""
    with open(os.path.join(model_dir, "install.yml")) as f:
        spec = yaml.safe_load(f)

    python_version = str(spec["python"])
    steps = []
    for entry in spec["commands"]:
        tool, name, *rest = entry
        if tool == "pip":
            is_vcs = name.startswith(("git+", "http://", "https://"))
            if is_vcs or not rest:
                package_spec, extra_args = name, rest
            else:
                package_spec, extra_args = f"{name}=={rest[0]}", rest[1:]
            steps.append(["pip", "install", package_spec, *extra_args])
        elif tool == "conda":
            rest = [x for x in rest if x not in ("-y", "--yes")]
            package_spec = f"{name}={rest[0]}" if rest else name
            cmd = ["conda", "install"]
            if len(rest) >= 2:
                cmd += ["-c", rest[1]]
            cmd += [package_spec, "-y"]
            steps.append(cmd)
        else:
            raise ValueError(f"install.yml: unknown command tool {tool!r}")
    return python_version, steps


def load_dockerfile(model_dir: str) -> tuple[str, list[str]]:
    """Parse Dockerfile into (python_version, [shell command, ...])."""
    with open(os.path.join(model_dir, "Dockerfile")) as f:
        text = f.read()

    from_line = next(
        (line for line in text.splitlines() if line.strip().startswith("FROM")), ""
    )
    match = DOCKERFILE_PY_RE.search(from_line)
    if not match:
        raise ValueError(
            f"Dockerfile: could not find a python version in {from_line!r}"
        )
    python_version = f"{match.group(1)}.{match.group(2)}"

    lines = text.replace("\\\n", " ").splitlines()
    commands = [
        line.strip()[len("RUN ") :].strip()
        for line in lines
        if line.strip().startswith("RUN ")
    ]
    return python_version, commands


def resolve_spec(model_dir: str) -> tuple[str, list]:
    """Dispatch to the install.yml or Dockerfile parser, whichever is present."""
    if os.path.exists(os.path.join(model_dir, "install.yml")):
        return load_install_yml(model_dir)
    if os.path.exists(os.path.join(model_dir, "Dockerfile")):
        return load_dockerfile(model_dir)
    raise FileNotFoundError(f"{model_dir}: no install.yml or Dockerfile found")


def run_step(prefix: str, step) -> subprocess.CompletedProcess:
    """Run one install step (argv list from install.yml, or a shell string
    from a Dockerfile RUN line) inside the env at `prefix`."""
    if isinstance(step, str):
        if CONDA_INSTALL_RE.match(step) and not CONDA_YES_RE.search(step):
            step = f"{step} -y"
        argv = ["conda", "run", "-p", prefix, "bash", "-c", step]
    else:
        argv = ["conda", "run", "-p", prefix, *step]
    return subprocess.run(
        argv, capture_output=True, text=True, check=False, env=ISOLATED_ENV
    )


def missing_requirements(prefix: str) -> list[str]:
    """Return `pip check` lines saying a requirement is not installed in the
    env (user site off). Version conflicts are not reported: some models pin
    conflicting versions on purpose."""
    result = subprocess.run(
        [os.path.join(prefix, "bin", "python"), "-m", "pip", "check"],
        capture_output=True,
        text=True,
        check=False,
        env=ISOLATED_ENV,
    )
    return [
        line.strip()
        for line in result.stdout.splitlines()
        if "which is not installed" in line
    ]


def build_env(identifier: str, models_dir: str, envs_dir: str, force: bool) -> bool:
    """Build the CPU conda env for one model. Returns True on full success."""
    model_dir = os.path.join(models_dir, identifier)
    prefix = os.path.join(envs_dir, f"{identifier}-cpu")

    if os.path.isdir(prefix):
        if not force:
            print(f"  {identifier}: already built, skipping")
            return True
        print(f"  {identifier}: --force, rebuilding")
        shutil.rmtree(prefix)

    try:
        python_version, steps = resolve_spec(model_dir)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"  {identifier}: could not resolve dependency spec - {exc}")
        return False

    print(f"  {identifier}: creating env (python={python_version})")
    result = subprocess.run(
        ["conda", "create", "-p", prefix, f"python={python_version}", "-y"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print(f"  {identifier}: conda create failed - {result.stderr.strip()}")
        return False

    for i, step in enumerate(steps, 1):
        shown = step if isinstance(step, str) else " ".join(step)
        print(f"  {identifier}: [{i}/{len(steps)}] {shown}")
        result = run_step(prefix, step)
        if result.returncode != 0:
            print(f"  {identifier}: step {i} failed - {result.stderr.strip()}")
            return False

    for line in missing_requirements(prefix):
        print(f"  {identifier}: WARNING {line}")

    print(f"  {identifier}: env ready at {prefix}")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--path-to-models",
        default=DEFAULT_MODELS_DIR,
        help=f"Where cloned models live (default: {DEFAULT_MODELS_DIR})",
    )
    parser.add_argument(
        "--path-to-envs",
        default=DEFAULT_ENVS_DIR,
        help=f"Where to build envs into (default: {DEFAULT_ENVS_DIR})",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="One or more model identifiers to build (default: every "
        "subdirectory of --path-to-models)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete and rebuild an env that already exists",
    )
    args = parser.parse_args()

    if shutil.which("conda") is None:
        raise SystemExit("conda not found on PATH")

    os.makedirs(args.path_to_envs, exist_ok=True)
    identifiers = find_models(args.path_to_models, args.models)
    print(f"Building {len(identifiers)} env(s) into {args.path_to_envs}")

    n_ok = 0
    for identifier in identifiers:
        if build_env(identifier, args.path_to_models, args.path_to_envs, args.force):
            n_ok += 1

    print(f"Done: {n_ok}/{len(identifiers)} environments ready")


if __name__ == "__main__":
    main()
