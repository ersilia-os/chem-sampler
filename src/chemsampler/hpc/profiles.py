"""Cluster profiles: how `chemsampler run --hpc NAME` reaches a Slurm cluster.

They live in `hpc.json` in chemsampler's folder (`$CHEMSAMPLER_HOME`, or
`~/.chemsampler`), a file of the user, never of a project or of chemsampler
itself, since it names hosts and nodes:

    {
      "mycluster": {
        "ssh": "login.example.org",
        "partition": "cpu",
        "nodes": ["node01", "node02"],
        "cpus": 4, "mem": "8G", "time": "01:00:00",
        "models": {"eos8nyg": {"mem": "16G"}}
      }
    }

Every key is optional. Leave `ssh` out when chemsampler runs where `sbatch` does.
"""

from dataclasses import dataclass, field
from pathlib import Path

from .. import store
from ..utils.jsonfile import is_list_of, is_object, is_text, read_object

HPC_FILE = "hpc.json"

#: What each job asks Slurm for; a profile's `models` can change them per model.
RESOURCES = ("cpus", "mem", "time")

KEYS = ("ssh", "partition", "nodes", *RESOURCES, "poll_seconds", "models")


@dataclass(frozen=True)
class Profile:
    """One cluster, as `hpc.json` describes it."""

    name: str
    ssh: str | None = None
    partition: str | None = None
    #: Passed as `--nodelist` with `--nodes=1`: a job may run on any of them.
    nodes: tuple[str, ...] = ()
    cpus: int = 1
    mem: str | None = None
    time: str | None = None
    #: How often the jobs are checked while a run waits for them.
    poll_seconds: float = 60.0
    #: Per-model resources: model id (or `mw`, `qed`) -> {resource: value}.
    models: dict = field(default_factory=dict)

    def resources(self, model_id: str) -> dict:
        """`cpus`, `mem` and `time` for one model: its own, else the profile's."""
        own = self.models.get(model_id, {})
        return {key: own.get(key, getattr(self, key)) for key in RESOURCES}


def profiles_path() -> Path:
    """The file with the profiles, in chemsampler's folder beside the store."""
    return store.home() / HPC_FILE


def load_profile(name: str) -> Profile:
    """
    Read one profile from `hpc.json` in chemsampler's folder.

    Raises
    ------
    ValueError
        If the file is missing or malformed, or has no profile of that name.
    """
    path = profiles_path()
    if not path.is_file():
        raise ValueError(
            f"no cluster profiles yet: write {path} (see docs/hpc.md in chemsampler)"
        )
    data = read_object(path, expected=" of profiles by name")
    if name not in data:
        known = ", ".join(sorted(data)) or "none"
        raise ValueError(f"{path} has no profile {name!r}; it has: {known}")
    return _parse(name, data[name], f"{path}: profile {name!r}")


def _parse(name: str, data, where: str) -> Profile:
    if not is_object(data):
        raise ValueError(f"{where} must be a JSON object")
    unknown = sorted(set(data) - set(KEYS))
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {unknown}; allowed: {list(KEYS)}")
    for key in ("ssh", "partition"):
        if key in data and not _filled(data[key]):
            raise ValueError(f"{where}: '{key}' must be text")
    if "nodes" in data and not (
        is_list_of(data["nodes"], str) and all(map(_filled, data["nodes"]))
    ):
        raise ValueError(f"{where}: 'nodes' must be a list of node names")
    if "poll_seconds" in data and not (
        isinstance(data["poll_seconds"], (int, float))
        and not isinstance(data["poll_seconds"], bool)
        and data["poll_seconds"] >= 0
    ):
        raise ValueError(f"{where}: 'poll_seconds' must be a number, 0 or more")
    _check_resources(data, where)
    models = data.get("models", {})
    if not is_object(models):
        raise ValueError(f"{where}: 'models' must map model ids to resources")
    for model_id, own in models.items():
        inner = f"{where}: models.{model_id}"
        if not is_object(own):
            raise ValueError(f"{inner} must be a JSON object")
        extra = sorted(set(own) - set(RESOURCES))
        if extra:
            raise ValueError(
                f"{inner}: unknown key(s) {extra}; allowed: {list(RESOURCES)}"
            )
        _check_resources(own, inner)
    return Profile(
        name=name,
        ssh=data.get("ssh"),
        partition=data.get("partition"),
        nodes=tuple(data.get("nodes", ())),
        cpus=data.get("cpus", 1),
        mem=data.get("mem"),
        time=data.get("time"),
        poll_seconds=float(data.get("poll_seconds", 60.0)),
        models=models,
    )


def _check_resources(data: dict, where: str) -> None:
    cpus = data.get("cpus", 1)
    if isinstance(cpus, bool) or not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"{where}: 'cpus' must be a whole number, 1 or more")
    for key in ("mem", "time"):
        if key in data and not _filled(data[key]):
            raise ValueError(f"{where}: '{key}' must be text, as sbatch takes it")


def _filled(value) -> bool:
    return is_text(value) and value.strip() != ""
