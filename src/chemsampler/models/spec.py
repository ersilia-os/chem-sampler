import math
from dataclasses import dataclass
from typing import Literal

Direction = Literal["higher", "lower"]


@dataclass(frozen=True)
class AnnotatorSpec:
    """
    An annotator paired with the cutoff a candidate must clear.

    Every annotator is treated uniformly: it has a cutoff and a direction, and
    nothing distinguishes one annotator from another beyond that. How annotators
    combine into a round's winner is decided by `hill_climb`'s `mode` argument,
    not by anything on the spec itself.

    Parameters
    ----------
    annotator_id : str
        Identifier used as the value's column name in round output (e.g. "qed",
        an Ersilia model id).
    annotator : object
        Any object exposing `.score(smiles_list: list[str]) -> dict[str, float]`.
    cutoff : float
        The value a candidate must clear to satisfy this annotator. Use
        `float("-inf")`/`float("inf")` to express "no real floor/ceiling".
    direction : {"higher", "lower"}
        "higher" is satisfied by `value >= cutoff`; "lower" by `value <= cutoff`.
    weight : float, optional
        Relative importance in `hill_climb`'s `mode="weighted"`; ignored by
        every other mode. Must be >= 0. By default 1.0.
    """

    annotator_id: str
    annotator: object
    cutoff: float
    direction: Direction
    weight: float = 1.0

    def __post_init__(self) -> None:
        if self.direction not in ("higher", "lower"):
            raise ValueError(
                f"direction must be 'higher' or 'lower', got {self.direction!r}"
            )
        if math.isnan(self.cutoff):
            raise ValueError(f"{self.annotator_id}: cutoff cannot be NaN")
        if math.isnan(self.weight):
            raise ValueError(f"{self.annotator_id}: weight cannot be NaN")
        if self.weight < 0:
            raise ValueError(
                f"{self.annotator_id}: weight must be >= 0, got {self.weight}"
            )
