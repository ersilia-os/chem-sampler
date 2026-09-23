from dataclasses import dataclass
from typing import Literal

Role = Literal["directing", "controlling"]


@dataclass(frozen=True)
class AnnotatorSpec:
    """
    An annotator paired with how it should be used in a round.

    A "directing" annotator is the optimization target: exactly one is required
    per `hill_climb` run, and its value becomes a round's `score`. A "controlling"
    annotator is a constraint: a candidate outside `[min, max]` is kept in the
    output but cannot win the round.

    Parameters
    ----------
    annotator_id : str
        Identifier used as the value's column name in round output (e.g. "qed",
        an Ersilia model id).
    annotator : object
        Any object exposing `.score(smiles_list: list[str]) -> dict[str, float]`.
    role : {"directing", "controlling"}
        Whether this annotator sets the optimization goal or a constraint.
    min : float, optional
        Lower bound for a controlling annotator. `None` means unbounded below.
    max : float, optional
        Upper bound for a controlling annotator. `None` means unbounded above.
    """

    annotator_id: str
    annotator: object
    role: Role
    min: float | None = None
    max: float | None = None

    def __post_init__(self) -> None:
        if self.role not in ("directing", "controlling"):
            raise ValueError(
                f"role must be 'directing' or 'controlling', got {self.role!r}"
            )
        if self.role == "directing" and (self.min is not None or self.max is not None):
            raise ValueError(
                f"{self.annotator_id}: a directing annotator cannot have min/max"
            )
        if self.role == "controlling" and self.min is None and self.max is None:
            raise ValueError(
                f"{self.annotator_id}: a controlling annotator needs at least one of min/max"
            )
