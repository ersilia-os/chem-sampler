"""Where a model's conda environment lives and whether its build finished."""

import json
from pathlib import Path

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


def read_marker(env: Path) -> dict | None:
    """What `write_marker` recorded, or None if the build did not finish."""
    try:
        return json.loads((Path(env) / ENV_MARKER).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None
