"""The workspace: one folder holding a run's inputs, models, environments and results."""

import json
from dataclasses import dataclass
from pathlib import Path

#: The marker file. Commands look for it in the working directory and its parents,
#: the way git looks for `.git`, and use it as their config file.
WORKSPACE_FILE = "chemsampler.json"

#: The keys of `WORKSPACE_FILE` and the folder each one defaults to.
FOLDER_KEYS = {
    "path_to_models": "input",
    "path_to_model_repositories": "models",
    "path_to_envs": "envs",
    "path_to_rounds": "rounds",
}
LOGS_FOLDER = "logs"

_GITIGNORE = "models/\nenvs/\nrounds/\nlogs/\n"

#: A first set of inputs, one model of each kind (same as `examples/inputs/csv`).
TEMPLATES = {
    "similarity.csv": "similarity\neos8nyg\n",
    "generators.csv": "generators\neos9taz\n",
    "annotators.csv": (
        "source,model,cutoff,direction,priority\nersilia,eos4zfy,0.0,higher,primary\n"
    ),
}


def _is_object(value) -> bool:
    return isinstance(value, dict)


@dataclass(frozen=True)
class Workspace:
    """The folders of a workspace, as its marker file says."""

    root: Path
    input: Path
    models: Path
    envs: Path
    rounds: Path
    logs: Path

    @classmethod
    def load(cls, root) -> "Workspace":
        """Read `root`'s marker file; folders it does not name get their defaults."""
        root = Path(root)
        path = root / WORKSPACE_FILE
        try:
            data = json.loads(path.read_text())
        except json.JSONDecodeError as e:
            raise ValueError(f"{path}: not valid JSON ({e})") from None
        if not _is_object(data):
            raise ValueError(f"{path}: expected a JSON object")

        def folder(key):
            value = Path(data.get(key, FOLDER_KEYS[key])).expanduser()
            return value if value.is_absolute() else root / value

        return cls(
            root=root,
            input=folder("path_to_models"),
            models=folder("path_to_model_repositories"),
            envs=folder("path_to_envs"),
            rounds=folder("path_to_rounds"),
            logs=root / LOGS_FOLDER,
        )


def find_workspace(start=None) -> Path | None:
    """The folder holding the marker file, looking from `start` (default: the cwd) upwards."""
    here = Path(start or Path.cwd()).resolve()
    for folder in (here, *here.parents):
        if (folder / WORKSPACE_FILE).is_file():
            return folder
    return None


def create_workspace(root) -> list[str]:
    """
    Create a workspace in `root`, adding whatever is missing and changing nothing else.

    Writes the marker file, the folders, a `.gitignore` for the bulky ones and a
    template of each input CSV. Anything that already exists, including inputs the
    user has edited, is left untouched, so running it again is safe.

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

    root.mkdir(parents=True, exist_ok=True)
    write_if_absent(root / WORKSPACE_FILE, json.dumps(FOLDER_KEYS, indent=2) + "\n")
    workspace = Workspace.load(root)
    for folder in (
        workspace.input,
        workspace.models,
        workspace.envs,
        workspace.rounds,
        workspace.logs,
    ):
        if not folder.exists():
            folder.mkdir(parents=True)
            created.append(f"{folder.relative_to(root)}/")
    write_if_absent(root / ".gitignore", _GITIGNORE)
    for name, text in TEMPLATES.items():
        write_if_absent(workspace.input / name, text)
    return created
