"""A project: one analysis's inputs and results, and where its models are.

The models themselves are not here. They are expensive to fetch and build, so by
default they live in the machine's store (see `chemsampler.store`), shared by
every project; the marker file records the folders. Switching project is `cd`.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from .utils.jsonfile import read_object

#: The marker file. Commands look for it in the working directory and its
#: parents, the way git looks for `.git`, and read it as their config file.
PROJECT_FILE = "chemsampler.json"

#: The folders a project owns, and the key that can move each one.
FOLDER_KEYS = {"path_to_models": "input", "path_to_rounds": "rounds"}

#: The folders of the models and their environments. `chemsampler init` always
#: writes them (the store's, unless told otherwise); without them, the store's.
STORE_KEYS = ("path_to_model_repositories", "path_to_envs")

_GITIGNORE = "rounds/\n"

#: A first set of inputs, one model of each kind (same as `examples/inputs/csv`).
TEMPLATES = {
    "similarity.csv": "similarity\neos8nyg\n",
    "generators.csv": "generators\neos9taz\n",
    "annotators.csv": (
        "source,model,cutoff,direction,priority\nersilia,eos4zfy,0.0,higher,primary\n"
    ),
}


@dataclass(frozen=True)
class Project:
    """The folders of a project, as its marker file says."""

    root: Path
    input: Path
    rounds: Path
    #: Only when the marker file names them; otherwise the store's are used.
    models: Path | None
    envs: Path | None

    @classmethod
    def load(cls, root) -> "Project":
        """Read `root`'s marker file; folders it does not name get their defaults."""
        root = Path(root)
        data = read_object(root / PROJECT_FILE)

        def folder(key, default=None):
            named = data.get(key, default)
            if named is None:
                return None
            value = Path(named).expanduser()
            return value if value.is_absolute() else root / value

        return cls(
            root=root,
            input=folder("path_to_models", FOLDER_KEYS["path_to_models"]),
            rounds=folder("path_to_rounds", FOLDER_KEYS["path_to_rounds"]),
            models=folder("path_to_model_repositories"),
            envs=folder("path_to_envs"),
        )


def find_project(start=None) -> Path | None:
    """The folder holding the marker file, looking from `start` (default: the cwd) upwards."""
    here = Path(start or Path.cwd()).resolve()
    for folder in (here, *here.parents):
        if (folder / PROJECT_FILE).is_file():
            return folder
    return None


def locate_project(path=None) -> Path:
    """
    The project that `path` (default: the working directory) is in.

    Raises
    ------
    ValueError
        If there is no `chemsampler.json` there or in any parent folder.
    """
    root = find_project(path)
    if root is not None:
        return root
    if path is None:
        raise ValueError(
            f"not inside a project (no {PROJECT_FILE} here or above); create one "
            "with `chemsampler init`, or give --path"
        )
    raise ValueError(
        f"{path} is not in a project (no {PROJECT_FILE} there or above); create "
        "one with `chemsampler init`"
    )


def create_project(root, models=None, envs=None) -> list[str]:
    """
    Create a project in `root`, adding whatever is missing and changing nothing else.

    Writes the marker file, the two folders it owns, a `.gitignore` for the
    results and a template of each input CSV. Anything that already exists,
    including inputs the user has edited, is left untouched, so running it again
    is safe.

    Parameters
    ----------
    root : str or Path
        The project folder.
    models, envs : str or Path, optional
        The folders of the model clones and their environments, recorded in a new
        marker file. Not given, the marker leaves them out and the store's are used.

    Returns
    -------
    list[str]
        What was created, as paths relative to `root` (empty if all was there).
    """
    root = Path(root)
    created = []

    def write_if_absent(path: Path, text: str):
        if not path.exists():
            path.write_text(text)
            created.append(str(path.relative_to(root)))

    marker = dict(FOLDER_KEYS)
    for key, folder in zip(STORE_KEYS, (models, envs)):
        if folder is not None:
            marker[key] = str(folder)
    root.mkdir(parents=True, exist_ok=True)
    write_if_absent(root / PROJECT_FILE, json.dumps(marker, indent=2) + "\n")
    project = Project.load(root)
    for folder in (project.input, project.rounds):
        if not folder.exists():
            folder.mkdir(parents=True)
            created.append(f"{folder.relative_to(root)}/")
    write_if_absent(root / ".gitignore", _GITIGNORE)
    for name, text in TEMPLATES.items():
        write_if_absent(project.input / name, text)
    return created
