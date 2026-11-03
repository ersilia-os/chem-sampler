"""Deleting models from the store, for `chemsampler remove`.

A model is a clone, an environment and an entry in fetch's record file. Nothing
here reads inside an environment, and the project's input lists are never touched.
"""

import shutil
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from .envs import RDKIT_ENV, env_path
from .fetch import _has_local_commits, _load_state, _local_changes, _run, _save_state
from .listing import present_models


@dataclass(frozen=True)
class Removal:
    """What `remove` would delete for one model."""

    model_id: str
    clone: Path | None
    env: Path | None
    #: Why it is kept (local work in the clone), or None.
    blocked: str | None = None

    @property
    def parts(self) -> str:
        """What goes, in words: "clone, environment"."""
        return ", ".join(
            name
            for name, there in (("clone", self.clone), ("environment", self.env))
            if there is not None
        )


@dataclass(frozen=True)
class RemovalPlan:
    models_dir: Path
    envs_dir: Path
    removals: tuple[Removal, ...]
    #: Ids that were asked for and are not in the store at all.
    absent: tuple[str, ...] = ()

    @property
    def to_remove(self) -> tuple[Removal, ...]:
        return tuple(r for r in self.removals if r.blocked is None)

    @property
    def blocked(self) -> tuple[Removal, ...]:
        return tuple(r for r in self.removals if r.blocked is not None)


def local_work(clone: Path) -> str | None:
    """
    Why a clone holds work that deleting would lose, or None. A folder that is
    no git clone (such as one `ersilia fetch` left) has nothing to compare.
    """
    if not (Path(clone) / ".git").is_dir():
        return None
    changed = _run(["git", "-C", clone, "status", "--porcelain"]).stdout
    if _local_changes(changed):
        return "local changes"
    if _has_local_commits(clone):
        return "local commits"
    return None


def _removal(model_id, models_dir: Path, envs_dir: Path, force: bool) -> Removal | None:
    """The removal of one model, or None if there is nothing of it on disk."""
    if model_id == RDKIT_ENV:  # an environment only: there is no clone to a built-in
        clone = None
    else:
        clone = models_dir / model_id
        clone = clone if clone.is_dir() else None
    env = env_path(envs_dir, model_id)
    env = env if env.is_dir() else None
    if clone is None and env is None:
        return None
    blocked = None if force or clone is None else local_work(clone)
    return Removal(model_id, clone, env, blocked)


def plan_removal(
    models_dir, envs_dir, ids: Sequence[str], force: bool = False
) -> RemovalPlan:
    """What removing these ids would delete: each one's clone and environment."""
    models_dir, envs_dir = Path(models_dir), Path(envs_dir)
    removals, absent = [], []
    for model_id in ids:
        found = _removal(model_id, models_dir, envs_dir, force)
        if found is None:
            absent.append(model_id)
        else:
            removals.append(found)
    return RemovalPlan(models_dir, envs_dir, tuple(removals), tuple(absent))


def _env_names(envs_dir: Path) -> set[str]:
    try:
        return {
            p.name
            for p in envs_dir.iterdir()
            if p.is_dir() and not p.name.startswith(".")
        }
    except (FileNotFoundError, NotADirectoryError):
        return set()


def plan_unused(
    models_dir,
    envs_dir,
    used: Sequence[str],
    rdkit_used: bool,
    force: bool = False,
) -> RemovalPlan:
    """
    What removing everything the project does not use would delete: every clone
    and environment whose model is not in `used`, orphan environments included,
    and the rdkit environment when no annotator needs it.
    """
    models_dir, envs_dir = Path(models_dir), Path(envs_dir)
    keep = set(used) | ({RDKIT_ENV} if rdkit_used else set())
    names = (present_models(models_dir) | _env_names(envs_dir)) - keep
    return plan_removal(models_dir, envs_dir, sorted(names), force)


def remove(
    plan: RemovalPlan, say: Callable[[str], None] = print
) -> list[tuple[str, str]]:
    """
    Delete what the plan says, one model at a time, skipping blocked ones.

    A failure on one model is returned as `(model_id, reason)` and does not stop
    the others. The record fetch keeps of a model goes with its clone.

    Returns
    -------
    list[tuple[str, str]]
        The models that could not be fully removed, with the reason.
    """
    failures, forget = [], []
    for item in plan.to_remove:
        say(f"{item.model_id}: removing ({item.parts})")
        try:
            if item.clone is not None:
                shutil.rmtree(item.clone)
                forget.append(item.model_id)  # its record goes with the clone
            if item.env is not None:
                shutil.rmtree(item.env)
        except OSError as e:
            failures.append((item.model_id, str(e)))
    state = _load_state(plan.models_dir)
    if any(model_id in state for model_id in forget):
        for model_id in forget:
            state.pop(model_id, None)
        _save_state(plan.models_dir, state)
    return failures
