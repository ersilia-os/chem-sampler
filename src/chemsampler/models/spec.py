import math
from dataclasses import dataclass
from typing import Literal

Direction = Literal["higher", "lower"]
AnnotatorDirection = Literal["higher", "lower", "range"]
Priority = Literal["primary", "secondary"]
Cutoff = float | tuple[float, float]


@dataclass(frozen=True)
class AnnotatorSpec:
    """
    An annotator paired with the cutoff a candidate must clear.

    A *primary* annotator is optimized: `hill_climb`'s `mode` decides how the
    primaries combine into a round's winner. A *secondary* annotator is a
    guard-rail: a candidate that does not satisfy it is never eligible, in any
    mode, and it takes no part in the optimization.

    Parameters
    ----------
    annotator_id : str
        Identifier used as the value's column name in round output (e.g. "qed",
        an Ersilia model id).
    annotator : object
        Any object exposing `.score(smiles_list: list[str]) -> dict[str, float]`.
    cutoff : float or tuple[float, float]
        The value a candidate must clear to satisfy this annotator. Use
        `float("-inf")`/`float("inf")` to express "no real floor/ceiling". A
        `(min, max)` pair, with `min < max`, when `direction="range"`.
    direction : {"higher", "lower", "range"}
        "higher" is satisfied by `value >= cutoff`; "lower" by `value <= cutoff`;
        "range" by `min <= value <= max`. A range is only allowed on a secondary
        annotator.
    weight : float, optional
        Relative importance in `hill_climb`'s `mode="weighted"`; ignored by
        every other mode and by secondary annotators. Must be >= 0. By default
        1.0.
    priority : {"primary", "secondary"}, optional
        "primary" annotators are optimized; "secondary" ones are guard-rails.
        By default "primary".
    """

    annotator_id: str
    annotator: object
    cutoff: Cutoff
    direction: AnnotatorDirection
    weight: float = 1.0
    priority: Priority = "primary"

    def __post_init__(self) -> None:
        if self.direction not in ("higher", "lower", "range"):
            raise ValueError(
                f"direction must be 'higher', 'lower' or 'range', got {self.direction!r}"
            )
        if self.priority not in ("primary", "secondary"):
            raise ValueError(
                f"{self.annotator_id}: priority must be 'primary' or 'secondary', "
                f"got {self.priority!r}"
            )
        if self.direction == "range":
            self._check_range()
        elif isinstance(self.cutoff, tuple):
            raise ValueError(
                f"{self.annotator_id}: cutoff {self.cutoff} is a range, so "
                f"direction must be 'range'; got {self.direction!r}"
            )
        elif math.isnan(self.cutoff):
            raise ValueError(f"{self.annotator_id}: cutoff cannot be NaN")
        if math.isnan(self.weight):
            raise ValueError(f"{self.annotator_id}: weight cannot be NaN")
        if self.weight < 0:
            raise ValueError(
                f"{self.annotator_id}: weight must be >= 0, got {self.weight}"
            )

    def _check_range(self) -> None:
        if self.priority != "secondary":
            raise ValueError(
                f"{self.annotator_id}: range is only allowed for secondary "
                "annotators (guard-rails); a primary annotator is optimized "
                "toward a one-sided cutoff"
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
