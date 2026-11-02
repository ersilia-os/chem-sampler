"""The store: the clones and conda environments every project on this machine shares.

Models are expensive — gigabytes of checkpoints and minutes of building — so they
live in one place instead of inside each project. `chemsampler fetch` fills it and
every project reads from it. It is chemsampler's one folder: the cluster profiles
(`hpc.json`) live there too.
"""

import os
from dataclasses import dataclass
from pathlib import Path

#: The environment variable that moves the folder, e.g. to a filesystem a
#: cluster's nodes share.
HOME_VARIABLE = "CHEMSAMPLER_HOME"

#: What fetch downloaded, recorded beside the clones it describes (never in the
#: store root: two projects may keep their clones in different folders).
STATE_FILE = ".chemsampler-fetch.json"


@dataclass(frozen=True)
class Store:
    """Where this machine keeps its model clones and their environments."""

    root: Path
    models: Path
    envs: Path

    @property
    def logs(self) -> Path:
        """Where `fetch` writes what it ran, per model."""
        return self.root / "logs"


def home() -> Path:
    """chemsampler's folder: `$CHEMSAMPLER_HOME`, or `~/.chemsampler` without it."""
    named = os.environ.get(HOME_VARIABLE)
    return Path(named or "~/.chemsampler").expanduser()


def resolve() -> Store:
    """
    The store this machine uses: `models/` and `envs/` inside `home()`.

    Raises
    ------
    ValueError
        If the folder holds a `config.json`, which older versions read to move the
        store and which would otherwise be silently ignored.
    """
    root = home()
    old = root / "config.json"
    if old.exists():
        raise ValueError(
            f"{old} is no longer read: chemsampler keeps everything in one folder. "
            f"Set {HOME_VARIABLE} to the folder it names (or move that folder "
            "here), then delete it."
        )
    return Store(root=root, models=root / "models", envs=root / "envs")


def state_path(models_dir) -> Path:
    """Where fetch records what it downloaded for the clones in `models_dir`."""
    return Path(models_dir) / STATE_FILE
