"""What is in the store, and what the Hub offers, for `chemsampler catalog`.

Everything here is cheap on purpose: it reads directory names and one marker file
per model, never an environment's contents. `pip list` or `du` over an
environment took tens of seconds on a shared filesystem.
"""

import dataclasses
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from .catalog import ROLES, Catalog
from .client import SIMILARITY_MODEL_IDS
from .envs import env_path, is_complete, read_marker
from .fetch import LIVE, release_problem


@dataclass(frozen=True)
class LocalModel:
    """A model in the store, or one a project needs that the store lacks."""

    model_id: str
    role: str  # "generators" | "similarity" | "annotators" | "unknown"
    title: str
    cloned: bool
    #: "ok" | "unfinished" | "missing" | "stale" (built from another commit than
    #: the clone is at) | "-" when no environments are used
    env: str
    used: bool  # named in one of the current project's input lists
    version: str | None = None  # what the inputs ask for
    at: str | None = None  # what the clone is at: a tag, "live" or a commit
    #: Why a used model is not ready, or None. Only used models have one.
    problem: str | None = None


@dataclass(frozen=True)
class HubEntry:
    """A model the Hub offers, and whether this machine already has it."""

    model_id: str
    role: str
    title: str
    in_store: bool
    supported: bool


def present_models(models_dir) -> set[str]:
    """The names of the clones in a folder: its sub-folders, hidden ones aside."""
    try:
        return {
            entry.name
            for entry in Path(models_dir).iterdir()
            if entry.is_dir() and not entry.name.startswith(".")
        }
    except (FileNotFoundError, NotADirectoryError):
        return set()


def _git(model_dir, *args) -> str | None:
    found = subprocess.run(
        ["git", "-C", str(model_dir), *args],
        capture_output=True,
        text=True,
        check=False,
    )
    return found.stdout.strip() if found.returncode == 0 else None


def _is_git(model_dir) -> bool:
    return (Path(model_dir) / ".git").is_dir()


def stale_commits(envs_dir, models_dir, model_id: str) -> tuple[str, str] | None:
    """
    The commit a complete environment was built from and the one the clone is at,
    if they differ; None if they agree or either is unknown (a model that is no
    git clone, an environment that records no commit).
    """
    clone = Path(models_dir) / model_id
    if not _is_git(clone):
        return None
    built = (read_marker(env_path(envs_dir, model_id)) or {}).get("commit")
    head = _git(clone, "rev-parse", "HEAD")
    if not built or not head or head.startswith(built):
        return None
    return built, head[: len(built)]


def env_state(envs_dir, model_id: str, models_dir=None) -> str:
    """`ok` if the environment was built completely, `unfinished` if a build was
    interrupted, `missing` if there is none. With `models_dir`, a complete one built
    from another commit than its clone is at is `stale`. Reads the marker only,
    never the environment."""
    env = env_path(envs_dir, model_id)
    if not env.exists():
        return "missing"
    if not is_complete(env):
        return "unfinished"
    if models_dir is not None and stale_commits(envs_dir, models_dir, model_id):
        return "stale"
    return "ok"


def version_at(model_dir, wanted: str | None = None) -> str | None:
    """
    What a clone is at: `wanted` if it is there, else "live" on the main branch,
    else the release tag it sits on, else its commit. None for a folder that is no
    git clone (a model `ersilia fetch` left), which has no version to compare.
    """
    if not _is_git(model_dir):
        return None
    if wanted and release_problem(Path(model_dir), wanted) is None:
        return wanted
    if release_problem(Path(model_dir), LIVE) is None:
        return LIVE
    return _git(model_dir, "describe", "--tags", "--exact-match", "HEAD") or _git(
        model_dir, "rev-parse", "--short", "HEAD"
    )


def _problem(model, wanted, at, models_dir, envs_dir) -> str | None:
    """Why a used model is not ready, in words, or None."""
    if not model.cloned:
        return "not cloned"
    reasons = []
    if wanted and at and release_problem(Path(models_dir) / model.model_id, wanted):
        reasons.append(f"cloned at {at}, the inputs ask for {wanted}")
    if model.env in ("missing", "unfinished"):
        reasons.append(f"environment {model.env}")
    elif model.env == "stale":
        built, head = stale_commits(envs_dir, models_dir, model.model_id)
        reasons.append(f"environment built from {built}, the clone is at {head}")
    return "; ".join(reasons) or None


def list_local(
    models_dir,
    envs_dir,
    *,
    used: Sequence[str] = (),
    versions: Mapping[str, str] | None = None,
    catalog: Catalog | None = None,
    role: str | None = None,
) -> list[LocalModel]:
    """
    The models in the store, with the state of each one's environment.

    Parameters
    ----------
    models_dir, envs_dir : path
        The folders with the clones and the environments. `envs_dir` may be None
        for models that bring their own interpreter (as `ersilia fetch` leaves
        them), which have no environment to report.
    used : Sequence[str]
        Models named in the current project's input lists. Those the store lacks
        are listed too, so a missing one is visible.
    versions : Mapping[str, str], optional
        The version the inputs ask for, per used model: a clone at another one is
        reported in its `problem`.
    catalog : Catalog, optional
        Gives each model its role and title; without it they are "unknown".
    role : str, optional
        Keep only this role.
    """
    known = catalog.by_id() if catalog is not None else {}
    present = present_models(models_dir)
    models = []
    versions = versions or {}
    for model_id in sorted(present | set(used)):
        entry = known.get(model_id)
        cloned = model_id in present
        wanted = versions.get(model_id)
        at = version_at(Path(models_dir) / model_id, wanted) if cloned else None
        model = LocalModel(
            model_id=model_id,
            role=entry.role if entry else "unknown",
            title=entry.title if entry else "",
            cloned=cloned,
            env="-" if envs_dir is None else env_state(envs_dir, model_id, models_dir),
            used=model_id in used,
            version=wanted if model_id in used else None,
            at=at,
        )
        if model.used:
            model = dataclasses.replace(
                model, problem=_problem(model, wanted, at, models_dir, envs_dir)
            )
        models.append(model)
    if role is not None:
        models = [m for m in models if m.role == role]
    order = (*ROLES, "unknown")
    return sorted(models, key=lambda m: (order.index(m.role), m.model_id))


def list_hub(catalog: Catalog, models_dir, role: str | None = None) -> list[HubEntry]:
    """The models the catalog lists, marked when they are in the store.

    Only the models in `SIMILARITY_MODEL_IDS` can be a similarity model, so the
    others are marked unsupported rather than left to look usable.
    """
    present = present_models(models_dir)
    entries = [
        HubEntry(
            model_id=e.identifier,
            role=e.role,
            title=e.title,
            in_store=e.identifier in present,
            supported=e.role != "similarity" or e.identifier in SIMILARITY_MODEL_IDS,
        )
        for e in catalog.entries
        if role is None or e.role == role
    ]
    return sorted(entries, key=lambda e: (ROLES.index(e.role), e.model_id))
