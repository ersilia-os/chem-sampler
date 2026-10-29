"""What is in the store, and what the Hub offers, for `chemsampler catalog`.

Everything here is cheap on purpose: it reads directory names and one marker file
per model, never an environment's contents. `pip list` or `du` over an
environment took tens of seconds on a shared filesystem.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .catalog import ROLES, Catalog
from .client import SIMILARITY_MODEL_IDS
from .envs import env_path, is_complete


@dataclass(frozen=True)
class LocalModel:
    """A model in the store, or one a project needs that the store lacks."""

    model_id: str
    role: str  # "generators" | "similarity" | "annotators" | "unknown"
    title: str
    cloned: bool
    env: str  # "ok" | "unfinished" | "missing" | "-" when no environments are used
    needed: bool  # asked for by the current project's inputs


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


def env_state(envs_dir, model_id: str) -> str:
    """`ok` if the environment was built completely, `unfinished` if a build was
    interrupted, `missing` if there is none. Reads one file, never the environment."""
    env = env_path(envs_dir, model_id)
    if not env.exists():
        return "missing"
    return "ok" if is_complete(env) else "unfinished"


def list_local(
    models_dir,
    envs_dir,
    *,
    needed: Sequence[str] = (),
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
    needed : Sequence[str]
        Models the current project's inputs ask for. Those the store lacks are
        listed too, so a missing one is visible.
    catalog : Catalog, optional
        Gives each model its role and title; without it they are "unknown".
    role : str, optional
        Keep only this role.
    """
    known = catalog.by_id() if catalog is not None else {}
    present = present_models(models_dir)
    models = []
    for model_id in sorted(present | set(needed)):
        entry = known.get(model_id)
        models.append(
            LocalModel(
                model_id=model_id,
                role=entry.role if entry else "unknown",
                title=entry.title if entry else "",
                cloned=model_id in present,
                env="-" if envs_dir is None else env_state(envs_dir, model_id),
                needed=model_id in needed,
            )
        )
    if role is not None:
        models = [m for m in models if m.role == role]
    order = (*ROLES, "unknown")
    return sorted(models, key=lambda m: (order.index(m.role), m.model_id))


def list_hub(catalog: Catalog, models_dir, role: str | None = None) -> list[HubEntry]:
    """The models the catalog lists, marked when they are in the store.

    Only `eos8nyg` can be a similarity model, so the others are marked
    unsupported rather than left to look usable.
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
