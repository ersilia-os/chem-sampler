"""The store: the clones and conda environments every project on this machine shares.

Models are expensive — gigabytes of checkpoints and minutes of building — so they
live in one place instead of inside each project. `chemsampler fetch` fills it and
every project reads from it.
"""

from dataclasses import dataclass
from pathlib import Path

from .utils.jsonfile import is_text, read_object

#: Where the store goes unless something says otherwise. Read through this module
#: (never copied into another one) so a test can point it somewhere harmless.
DEFAULT_STORE = Path.home() / ".chemsampler"

#: The folder each setting names, inside the store.
STORE_FOLDERS = {"path_to_model_repositories": "models", "path_to_envs": "envs"}

#: What `<default store>/config.json` may hold: the store's own location, and a
#: folder-by-folder override of it.
CONFIG_KEYS = ("store", *STORE_FOLDERS)

CONFIG_FILE = "config.json"

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


def config_path() -> Path:
    """The machine-wide settings file, always next to the default store."""
    return DEFAULT_STORE / CONFIG_FILE


def read_config() -> dict:
    """
    Read the machine-wide settings, which need not exist.

    Returns
    -------
    dict
        The settings, or an empty dict when there is no file.

    Raises
    ------
    ValueError
        If the file is not a JSON object, has an unknown key, or a value that is
        not a path.
    """
    path = config_path()
    data = read_object(
        path, missing_ok=True, expected=f" with keys {list(CONFIG_KEYS)}"
    )
    unknown = sorted(set(data) - set(CONFIG_KEYS))
    if unknown:
        raise ValueError(
            f"{path}: unknown key(s) {unknown}; allowed: {list(CONFIG_KEYS)}"
        )
    for key, value in data.items():
        if not is_text(value):
            raise ValueError(f"{path}: '{key}' must be a path (text)")
    return data


def resolve() -> Store:
    """
    The store this machine uses.

    `<default store>/config.json` can move the whole store with `store`, or a
    single folder with `path_to_model_repositories` or `path_to_envs` — which is
    how a cluster puts the models on a shared filesystem.
    """
    data = read_config()
    root = Path(data["store"]).expanduser() if "store" in data else DEFAULT_STORE

    def folder(key: str) -> Path:
        if key in data:
            named = Path(data[key]).expanduser()
            return named if named.is_absolute() else root / named
        return root / STORE_FOLDERS[key]

    return Store(
        root=root,
        models=folder("path_to_model_repositories"),
        envs=folder("path_to_envs"),
    )


def state_path(models_dir) -> Path:
    """Where fetch records what it downloaded for the clones in `models_dir`."""
    return Path(models_dir) / STATE_FILE
