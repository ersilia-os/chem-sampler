"""Parsing and validation of the annotator inputs, shared by the CSV and JSON files."""

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..models.annotator import RDKIT_ANNOTATORS
from ..models.spec import Cutoff

SOURCES = ("ersilia", "rdkit")
DIRECTIONS = ("higher", "lower", "range")
PRIORITIES = ("primary", "secondary")
REQUIRED_FIELDS = ("source", "model", "cutoff", "direction", "priority")
OPTIONAL_FIELDS = ("weight", "column")

#: Two numbers joined by a hyphen, so "-5--1" is a range and "-5" is not.
_RANGE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*-\s*(-?\d+(?:\.\d+)?)\s*$")


@dataclass(frozen=True)
class AnnotatorRow:
    """One validated annotator entry, before its annotator object is built."""

    source: str
    model: str
    cutoff: Cutoff
    direction: str
    priority: str
    weight: float | None
    column: str | None


def check_header(
    path, header: Sequence[str] | None, required: Sequence[str], optional=()
) -> None:
    """Fail unless a CSV header holds every required column and nothing unknown."""
    columns = [c.strip() for c in header or []]
    missing = [c for c in required if c not in columns]
    unknown = [c for c in columns if c not in required and c not in optional]
    duplicated = sorted({c for c in columns if columns.count(c) > 1})
    if missing or unknown or duplicated:
        allowed = list(required) + list(optional)
        raise ValueError(
            f"{path}: bad header {columns}; expected {allowed}"
            + (f" (missing: {missing})" if missing else "")
            + (f" (unknown: {unknown})" if unknown else "")
            + (f" (repeated: {duplicated})" if duplicated else "")
        )


def parse_cutoff(where: str, direction: str, raw) -> Cutoff:
    """
    Turn a cutoff as written in a file into a number or a `(min, max)` pair.

    Parameters
    ----------
    where : str
        Prefix for error messages (the file, row and model).
    direction : {"higher", "lower", "range"}
        The annotator's direction, which decides whether a range is expected.
    raw : int, float or str
        A JSON number, or text: a number ("0.3", "1e-3") or "min-max"
        ("250-450", "-5--1").

    Raises
    ------
    ValueError
        If the cutoff is empty or unreadable, is a range but the direction is not
        "range", is a single number but the direction is "range", or a range's min
        is not lower than its max.
    """
    if (
        raw is None
        or isinstance(raw, bool)
        or (isinstance(raw, str) and not raw.strip())
    ):
        raise ValueError(f"{where}: cutoff is required")

    if isinstance(raw, (int, float)):
        cutoff: Cutoff = float(raw)
    else:
        text = str(raw).strip()
        try:
            cutoff = float(text)
        except ValueError:
            match = _RANGE.match(text)
            if match is None:
                raise ValueError(
                    f"{where}: invalid cutoff {raw!r}; write a number, or "
                    "'min-max' for a range"
                ) from None
            cutoff = (float(match[1]), float(match[2]))

    is_range = isinstance(cutoff, tuple)
    if direction == "range" and not is_range:
        raise ValueError(
            f"{where}: direction is 'range', so cutoff must be written as "
            f"'min-max' (e.g. '250-450'); got {raw!r}. Use direction 'higher' "
            "or 'lower' for a single threshold."
        )
    if is_range and direction != "range":
        raise ValueError(
            f"{where}: cutoff {raw!r} is a range, so direction must be "
            f"'range'; got {direction!r}."
        )
    if is_range and not cutoff[0] < cutoff[1]:
        raise ValueError(f"{where}: in the range {raw!r}, min must be lower than max")
    return cutoff


def _clean(value):
    """Strip text and treat an empty cell as missing."""
    if isinstance(value, str):
        value = value.strip()
        return value or None
    return value


def _parse_weight(where: str, value) -> float | None:
    if value is None:
        return None
    try:
        if isinstance(value, bool):
            raise TypeError("a boolean is not a weight")
        weight = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{where}: invalid weight {value!r}") from None
    if weight < 0:
        raise ValueError(f"{where}: weight must be >= 0, got {weight}")
    return weight


def parse_annotator(row: Mapping, where: str) -> AnnotatorRow:
    """
    Validate one annotator entry, from a CSV row or a JSON object.

    Parameters
    ----------
    row : Mapping
        Fields `source`, `model`, `cutoff`, `direction` and `priority`, and
        optionally `weight` and `column`. Empty cells count as missing.
    where : str
        Prefix for error messages, e.g. "annotators.csv row 3".

    Returns
    -------
    AnnotatorRow
        The validated entry; `weight` is `None` unless it was given.

    Raises
    ------
    ValueError
        On an unknown or missing field, a source/model/direction/priority that is
        not allowed, a bad cutoff (see `parse_cutoff`), a range on a primary
        annotator, a weight on a secondary one, or a column on an rdkit one.
    """
    allowed = REQUIRED_FIELDS + OPTIONAL_FIELDS
    unknown = sorted(str(k) for k in row if k not in allowed)
    if unknown:
        raise ValueError(
            f"{where}: unknown field(s) {unknown}; allowed: {list(allowed)}"
        )

    values = {name: _clean(row.get(name)) for name in allowed}
    missing = [name for name in REQUIRED_FIELDS if values[name] is None]
    if missing:
        raise ValueError(f"{where}: missing {missing}")

    source, model = values["source"], str(values["model"])
    direction, priority = values["direction"], values["priority"]
    where = f"{where} ({model})"

    if source not in SOURCES:
        raise ValueError(
            f"{where}: source must be one of {list(SOURCES)}, got {source!r}"
        )
    if source == "rdkit":
        if model not in RDKIT_ANNOTATORS:
            raise ValueError(
                f"{where}: unknown rdkit model; choose from {sorted(RDKIT_ANNOTATORS)}"
            )
        if values["column"] is not None:
            raise ValueError(f"{where}: column only applies to ersilia models")
    if direction not in DIRECTIONS:
        raise ValueError(
            f"{where}: direction must be one of {list(DIRECTIONS)}, got {direction!r}"
        )
    if priority not in PRIORITIES:
        raise ValueError(
            f"{where}: priority must be one of {list(PRIORITIES)}, got {priority!r}"
        )

    cutoff = parse_cutoff(where, direction, values["cutoff"])
    if direction == "range" and priority != "secondary":
        raise ValueError(
            f"{where}: range is only allowed for secondary annotators "
            "(guard-rails); a primary annotator is optimized toward a "
            "one-sided cutoff."
        )

    weight = _parse_weight(where, values["weight"])
    if priority == "secondary" and weight is not None:
        raise ValueError(
            f"{where}: weight only applies to primary annotators; a secondary "
            "annotator is a guard-rail and is not weighted."
        )

    return AnnotatorRow(
        source=source,
        model=model,
        cutoff=cutoff,
        direction=direction,
        priority=priority,
        weight=weight,
        column=values["column"],
    )
