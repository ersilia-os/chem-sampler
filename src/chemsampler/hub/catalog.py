"""The Hub labels chemsampler relies on, kept as a catalog file shipped with the package.

Runs check every model against `catalog.csv`, which is instant and works offline.
`chemsampler setup` rebuilds that file from the Hub with `ersilia catalog`; nothing
else queries the Hub.
"""

import csv
import datetime
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..utils.logging import logger
from .client import SIMILARITY_MODEL_IDS, not_allowed_similarity_message

#: The catalog file, in the repository and shipped with the package.
CATALOG_PATH = Path(__file__).with_name("catalog.csv")

ROLES = ("generators", "similarity", "annotators")
_COLUMNS = ["role", "identifier", "title", "task", "subtask"]
_HEADER = re.compile(
    r"^# chemsampler catalog, built (\d{4}-\d{2}-\d{2}), ersilia (\S+)$"
)
_AS = {
    "generators": "a generator",
    "similarity": "a similarity model",
    "annotators": "an annotator",
}
_SUBTASK_ROLE = {"Generation": "generators", "Similarity search": "similarity"}

_TIMEOUT_SECONDS = 60

#: What `ersilia catalog --card` prints for an id the Hub does not know. It
#: exits 0 and writes no file, exactly like an unreachable Hub does, so this text
#: is the only way to tell the two apart.
_UNKNOWN_ID_TEXT = "No metadata found"


@dataclass(frozen=True)
class CatalogEntry:
    """One Hub model and the role its labels give it."""

    role: str
    identifier: str
    title: str
    task: str
    subtask: str


@dataclass(frozen=True)
class Catalog:
    """The catalog file: when and with which ersilia it was built, and its models."""

    built: str
    ersilia_version: str
    entries: tuple[CatalogEntry, ...]

    def by_id(self) -> dict[str, CatalogEntry]:
        return {entry.identifier: entry for entry in self.entries}

    def count(self, role: str) -> int:
        return sum(entry.role == role for entry in self.entries)


@dataclass(frozen=True)
class CatalogDiff:
    """Models that appeared in, or left, the Hub lists since the file was built."""

    added: tuple[CatalogEntry, ...]
    removed: tuple[CatalogEntry, ...]

    @property
    def empty(self) -> bool:
        return not self.added and not self.removed


@dataclass(frozen=True)
class CatalogUpdate:
    """What `update_catalog` found and did."""

    catalog: Catalog
    previous: Catalog | None
    diff: CatalogDiff
    written: bool


# --- the file -----------------------------------------------------------------------


def read_catalog(path=None) -> Catalog:
    """
    Read the catalog file.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the catalog shipped with the package.

    Raises
    ------
    ValueError
        If the file is missing or is not a valid catalog.
    """
    path = Path(path or CATALOG_PATH)
    rebuild = "run `chemsampler setup` to build it"
    try:
        lines = path.read_text().splitlines()
    except FileNotFoundError:
        raise ValueError(
            f"the model catalog {path} does not exist; {rebuild}"
        ) from None
    header = _HEADER.match(lines[0]) if lines else None
    if header is None:
        raise ValueError(f"{path} is not a model catalog (bad first line); {rebuild}")
    reader = csv.DictReader(lines[1:])
    if reader.fieldnames != _COLUMNS:
        raise ValueError(f"{path}: expected columns {_COLUMNS}; {rebuild}")
    entries = []
    for number, row in enumerate(reader, start=1):
        if row["role"] not in ROLES or not row["identifier"]:
            raise ValueError(f"{path} row {number} is not valid; {rebuild}")
        entries.append(CatalogEntry(**row))
    return Catalog(header[1], header[2], tuple(entries))


def _today() -> str:
    """Today's date in UTC, so the same day is recorded wherever the file is rebuilt."""
    return datetime.datetime.now(tz=datetime.timezone.utc).date().isoformat()


def write_catalog(entries: Sequence[CatalogEntry], path=None) -> None:
    """
    Write the catalog file, dated today, replacing any previous one.

    Models are sorted by role and id so a rebuild changes only what changed.
    The file is written whole or not at all.
    """
    path = Path(path or CATALOG_PATH)
    try:
        version = importlib.metadata.version("ersilia")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    ordered = sorted(entries, key=lambda e: (ROLES.index(e.role), e.identifier))
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "w", newline="") as f:
        f.write(f"# chemsampler catalog, built {_today()}, ersilia {version}\n")
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(_COLUMNS)
        for e in ordered:
            writer.writerow([e.role, e.identifier, e.title, e.task, e.subtask])
    os.replace(temporary, path)


def diff_catalogs(old: Catalog | None, new: Sequence[CatalogEntry]) -> CatalogDiff:
    """Models in `new` but not in `old` (added) and in `old` but not in `new` (removed)."""
    old_entries = old.entries if old is not None else ()
    old_keys = {(e.role, e.identifier) for e in old_entries}
    new_keys = {(e.role, e.identifier) for e in new}
    return CatalogDiff(
        added=tuple(e for e in new if (e.role, e.identifier) not in old_keys),
        removed=tuple(e for e in old_entries if (e.role, e.identifier) not in new_keys),
    )


# --- the Hub (only `chemsampler setup` uses this) ------------------------------------


class _NoOutput(Exception):
    """`ersilia catalog` exited cleanly but wrote no file; carries what it printed."""


def _ersilia_executable() -> str:
    sibling = Path(sys.executable).with_name("ersilia")
    if sibling.exists():
        return str(sibling)
    found = shutil.which("ersilia")
    if found:
        return found
    raise RuntimeError(
        "the ersilia command was not found next to this Python or on PATH"
    )


def _unreachable(reason: str) -> RuntimeError:
    return RuntimeError(
        f"could not query the Ersilia Model Hub ({reason}); check your network connection"
    )


def _run_ersilia(*catalog_args: str):
    """Run `ersilia catalog --hub <args>` and return the JSON it writes."""
    with tempfile.TemporaryDirectory(prefix="chemsampler-") as tmp_dir:
        output = os.path.join(tmp_dir, "catalog.json")
        command = [
            _ersilia_executable(),
            "catalog",
            "--hub",
            *catalog_args,
            "--output",
            output,
        ]
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=_TIMEOUT_SECONDS,
                check=False,
            )
        except subprocess.TimeoutExpired:
            raise _unreachable(f"timed out after {_TIMEOUT_SECONDS} s") from None
        if result.returncode != 0:
            tail = (result.stderr or result.stdout).strip()[-300:]
            raise _unreachable(f"ersilia exited with {result.returncode}: {tail}")
        if not os.path.exists(output):
            raise _NoOutput((result.stdout or "") + (result.stderr or ""))
        try:
            with open(output) as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            raise _unreachable(f"unreadable ersilia output: {e}") from None


def _list_models(task: str) -> list[dict]:
    """The Hub models listed under one task (`sampling`, `annotation`)."""
    try:
        return _run_ersilia("--more", "--task", task)
    except _NoOutput as e:
        raise _unreachable(f"no catalog written: {str(e).strip()[-300:]}") from None


def _card(model_id: str) -> dict | None:
    """The model's full card, or None if the Hub does not know the id."""
    try:
        return _run_ersilia("--card", model_id)
    except _NoOutput as e:
        if _UNKNOWN_ID_TEXT in str(e):
            return None
        raise _unreachable(f"no card written: {str(e).strip()[-300:]}") from None


def build_catalog() -> list[CatalogEntry]:
    """
    Build the catalog from the Hub, which takes about 20 s.

    Two list calls give the Sampling and Annotation models, and one card per
    Sampling model gives the Subtask that separates generators from similarity
    search. Sampling models with any other Subtask are left out: no role uses them.

    Raises
    ------
    RuntimeError
        If the Hub cannot be queried, or lists a model it has no card for.
    """
    entries = []
    sampling = _list_models("sampling")
    logger.info(f"Reading the cards of {len(sampling)} Sampling models")
    for row in sampling:
        card = _card(row["Identifier"])
        if card is None:
            raise _unreachable(
                f"the Hub lists {row['Identifier']} but has no card for it"
            )
        role = _SUBTASK_ROLE.get(card.get("Subtask"))
        if role is not None:
            entries.append(
                CatalogEntry(
                    role, row["Identifier"], row["Title"], row["Task"], card["Subtask"]
                )
            )
    for row in _list_models("annotation"):
        entries.append(
            CatalogEntry("annotators", row["Identifier"], row["Title"], row["Task"], "")
        )
    return entries


def update_catalog(path=None) -> CatalogUpdate:
    """
    Rebuild the catalog from the Hub, rewriting the file only if something changed.

    A file that is missing or unreadable is rebuilt. If the Hub cannot be reached,
    the exception propagates and the existing file is left exactly as it was.

    Returns
    -------
    CatalogUpdate
        The catalog now on disk, the previous one (None if there was none), the
        models added and removed, and whether the file was written.
    """
    path = Path(path or CATALOG_PATH)
    try:
        previous = read_catalog(path)
    except ValueError:
        previous = None
    entries = build_catalog()
    diff = diff_catalogs(previous, entries)
    if previous is not None and diff.empty:
        return CatalogUpdate(previous, previous, diff, written=False)
    write_catalog(entries, path)
    return CatalogUpdate(read_catalog(path), previous, diff, written=True)


# --- validating models against the file ---------------------------------------------


def validate_models(
    *,
    generators: Sequence[str] = (),
    similarity: Sequence[str] = (),
    annotators: Sequence[str] = (),
    catalog: Catalog | None = None,
) -> None:
    """
    Check that every model carries the label its role needs, per the catalog file.

    Generators must be labelled Sampling / Generation, similarity models
    Sampling / Similarity search (and be one chemsampler supports), annotators
    the Annotation task. The catalog holds only models the Hub lists, which are
    the Ready ones.

    Parameters
    ----------
    generators : Sequence[str]
        Hub ids to be used as generators.
    similarity : Sequence[str]
        Hub ids to be used as similarity-search models.
    annotators : Sequence[str]
        Hub ids of annotators. In-house (rdkit) annotators are not Hub models
        and must not be passed.
    catalog : Catalog, optional
        Defaults to the catalog file shipped with the package.

    Raises
    ------
    ValueError
        If the catalog file is missing or invalid, or any model is not in it or
        has the wrong label. Every problem found is reported in the same error.
    """
    if not (generators or similarity or annotators):
        return
    catalog = catalog or read_catalog()
    logger.info(f"Checking models against the catalog built on {catalog.built}")
    known = catalog.by_id()

    def problem(model_id: str, role: str) -> str | None:
        entry = known.get(model_id)
        if entry is None:
            return (
                f"{model_id} ({role}): not in the catalog built on {catalog.built}; "
                "it is misspelled, not Ready, or newer than the catalog. Run "
                "`chemsampler setup` to refresh the catalog from the Hub."
            )
        supported = entry.role != "similarity" or model_id in SIMILARITY_MODEL_IDS
        if entry.role == role:
            return None if supported else not_allowed_similarity_message(model_id)
        label = f"{entry.task} / {entry.subtask}" if entry.subtask else entry.task
        message = f"{model_id} ({role}): is labelled {label}, so it cannot be used as {_AS[role]}"
        if supported:
            return f"{message}; put it in {entry.role}."
        return (
            f"{message}, and it is not meant to be used in chemsampler as a "
            f"similarity model either (choose from: "
            f"{', '.join(sorted(SIMILARITY_MODEL_IDS))})."
        )

    problems = [
        found
        for role, ids in (
            ("generators", generators),
            ("similarity", similarity),
            ("annotators", annotators),
        )
        for model_id in ids
        if (found := problem(model_id, role))
    ]
    if problems:
        raise ValueError("these models cannot be used:\n  - " + "\n  - ".join(problems))
