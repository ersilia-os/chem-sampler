import math
from dataclasses import dataclass
from typing import Literal

Direction = Literal["higher", "lower"]
AnnotatorDirection = Literal["higher", "lower", "range"]
Priority = Literal["primary", "guard-rail"]
Cutoff = float | tuple[float, float]


@dataclass(frozen=True)
class AnnotatorSpec:
    """
    An annotator paired with the cutoff a candidate must clear.

    A run has exactly one *primary* annotator, the one it optimizes. A
    *guard-rail* is a hard constraint: a candidate that does not satisfy it is
    never eligible, and it takes no part in the optimization.

    Parameters
    ----------
    annotator_id : str
        Identifier used as the value's column name in round output (e.g. "qed",
        an Ersilia model id).
    annotator : object
        Any object exposing `.score(smiles_list: list[str]) -> dict[str, float]`.
    cutoff : float or tuple[float, float]
        The value a candidate must clear to satisfy this annotator. A
        `(min, max)` pair, with `min < max`, when `direction="range"`. NaN, on a
        primary only, means no cutoff: the primary is simply optimized.
    direction : {"higher", "lower", "range"}
        "higher" is satisfied by `value >= cutoff`; "lower" by `value <= cutoff`;
        "range" by `min <= value <= max`. A range is only allowed on a guard-rail.
    priority : {"primary", "guard-rail"}, optional
        "primary" annotators are optimized; "guard-rail" ones are constraints.
        By default "primary".
    """

    annotator_id: str
    annotator: object
    cutoff: Cutoff
    direction: AnnotatorDirection
    priority: Priority = "primary"

    def __post_init__(self) -> None:
        if self.direction not in ("higher", "lower", "range"):
            raise ValueError(
                f"direction must be 'higher', 'lower' or 'range', got {self.direction!r}"
            )
        if self.priority not in ("primary", "guard-rail"):
            raise ValueError(
                f"{self.annotator_id}: priority must be 'primary' or 'guard-rail', "
                f"got {self.priority!r}"
            )
        if self.direction == "range":
            self._check_range()
        elif isinstance(self.cutoff, tuple):
            raise ValueError(
                f"{self.annotator_id}: cutoff {self.cutoff} is a range, so "
                f"direction must be 'range'; got {self.direction!r}"
            )
        elif math.isnan(self.cutoff) and self.priority == "guard-rail":
            raise ValueError(
                f"{self.annotator_id}: a guard-rail needs a cutoff (NaN, no "
                "cutoff, is only for the primary)"
            )

    def _check_range(self) -> None:
        if self.priority != "guard-rail":
            raise ValueError(
                f"{self.annotator_id}: range is only allowed for guard-rails; a "
                "primary annotator is optimized toward a one-sided cutoff"
            )
        if not (isinstance(self.cutoff, tuple) and len(self.cutoff) == 2):
            raise ValueError(
                f"{self.annotator_id}: direction 'range' needs a (min, max) "
                f"cutoff, got {self.cutoff!r}"
            )
        low, high = self.cutoff
        if not (math.isfinite(low) and math.isfinite(high) and low < high):
            raise ValueError(
                f"{self.annotator_id}: a range needs finite bounds with "
                f"min < max, got {self.cutoff}"
            )


def has_cutoff(spec: AnnotatorSpec) -> bool:
    """Whether the annotator has a cutoff; a primary without one (NaN) is optimized."""
    return not (isinstance(spec.cutoff, float) and math.isnan(spec.cutoff))
