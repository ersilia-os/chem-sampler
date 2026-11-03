"""Where a model's conda environment lives and whether its build finished."""

import json
import os
from pathlib import Path

#: The environment the in-house RDKit annotators (`mw`, `qed`, `tanimoto`) run
#: in, here and on a cluster: no model, only RDKit. `chemsampler fetch rdkit`
#: builds it.
RDKIT_ENV = "rdkit"

#: Written inside an environment when `chemsampler fetch` has built it completely.
#: An environment without it was interrupted or failed and is not to be trusted.
ENV_MARKER = ".chemsampler-env.json"


def env_path(envs_dir, model_id: str) -> Path:
    """The environment folder of one model."""
    return Path(envs_dir) / model_id


def is_complete(env: Path) -> bool:
    """True if `fetch` finished building the environment."""
    return (Path(env) / ENV_MARKER).is_file()


def write_marker(env: Path, **info) -> None:
    """Record that the environment was built completely, with what it was built from."""
    (Path(env) / ENV_MARKER).write_text(json.dumps(info, indent=2) + "\n")


def update_marker(env: Path, **info) -> None:
    """Change what the marker of a complete environment records, keeping the rest."""
    write_marker(env, **{**(read_marker(env) or {}), **info})


def read_marker(env: Path) -> dict | None:
    """What `write_marker` recorded, or None if the build did not finish."""
    try:
        return json.loads((Path(env) / ENV_MARKER).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def process_env(env_dir: Path) -> dict[str, str]:
    """
    The environment variables of a process run as if `env_dir` were active.

    A model cloned and built by `chemsampler fetch` has a plain `run.sh`
    (`python main.py ...`): the env's `bin` goes first on PATH and its own
    libraries first on the loader path (otherwise the system's older libstdc++
    is loaded), the user site is off, and nothing of chemsampler's own
    PYTHONPATH is kept. A cluster job sets the same (`chemsampler.hpc.jobs`).
    """
    env = {**os.environ}
    env.pop("PYTHONPATH", None)
    env["PATH"] = os.pathsep.join(filter(None, [str(env_dir / "bin"), env.get("PATH")]))
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        filter(None, [str(env_dir / "lib"), env.get("LD_LIBRARY_PATH")])
    )
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    return env
