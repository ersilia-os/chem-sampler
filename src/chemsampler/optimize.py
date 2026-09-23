import math

import pandas as pd

from .models.spec import AnnotatorSpec
from .utils.logging import logger
from .utils.similarity import tanimoto_to_seed

_RESERVED_COLUMNS = {"smiles", "source", "passes_constraints", "tanimoto_to_seed"}


def hill_climb(
    generator,
    annotators: list[AnnotatorSpec],
    *,
    seed_smiles: str | None = None,
    n_rounds: int = 5,
    tolerance: float = 0.0,
    tanimoto_cutoff: float | None = None,
) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    """
    Iteratively improve a seed molecule's score by alternating generation and scoring.

    Each round, `generator` proposes candidates from the current-best molecule and
    every annotator in `annotators` scores them. Among candidates that satisfy all
    controlling annotators (and the Tanimoto cutoff, if any), the best directing
    score becomes the seed for the next round. Stops after `n_rounds`, or as soon
    as a round produces no eligible candidate, or as soon as a round fails to beat
    the current best by more than `tolerance`.

    A candidate that fails a constraint is kept in the round's output (with
    `passes_constraints=False`) rather than dropped — it just cannot win the round.

    Parameters
    ----------
    generator : object
        Any object exposing `.generate_by_model(seed_smiles: str | None) ->
        dict[str, list[str]]` (e.g. `HubGenerator`, `ChemblSampler`, or a
        `GeneratorPool`). If `seed_smiles` is `None`, a generator that structurally
        requires a seed should raise `SeedRequiredError`; `GeneratorPool` already
        turns that into a warning and skips the generator.
    annotators : list[AnnotatorSpec]
        Exactly one must have `role="directing"` (the optimization target); any
        number may have `role="controlling"` (constraints).
    seed_smiles : str, optional
        SMILES of the starting molecule. If `None`, round 0 is skipped and the
        first round's winner is unconditionally the new best.
    n_rounds : int, optional
        Maximum number of rounds to run, by default 5.
    tolerance : float, optional
        A round only counts as an improvement if it beats the current best score
        by more than this margin, by default 0.0 (any improvement counts).
    tanimoto_cutoff : float, optional
        Minimum Tanimoto similarity to `seed_smiles` for a candidate to satisfy
        constraints. Requires `seed_smiles` to be given.

    Returns
    -------
    summary : pandas.DataFrame
        One row per round attempted, with columns `round`, `smiles` (that round's
        winner), `score` (its directing-annotator value), `is_new_best`.
    candidates_by_round : dict[int, pandas.DataFrame]
        Every candidate considered in each round (including round 0, the seed, if
        a seed was given), with columns `smiles`, `source` (contributing
        generator id(s), joined with `|` if more than one produced the same
        molecule), one column per `annotators` entry's `annotator_id`,
        `passes_constraints`, and `tanimoto_to_seed` (only present if a seed was
        given; always computed against the original seed, never the rolling
        best-so-far molecule).

    Raises
    ------
    ValueError
        If `annotators` doesn't have exactly one directing spec, has duplicate or
        reserved `annotator_id` values, `tanimoto_cutoff` is given without a seed,
        or the seed cannot be scored by the directing annotator.
    """
    directing = _validate_annotators(annotators)
    if tanimoto_cutoff is not None and seed_smiles is None:
        raise ValueError("tanimoto_cutoff requires seed_smiles")

    best_smiles = seed_smiles
    best_score = float("-inf")
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}

    if seed_smiles is not None:
        seed_scores = _score_candidates([seed_smiles], annotators)
        seed_table = _round_table(
            {seed_smiles: "seed"},
            seed_scores,
            annotators,
            {seed_smiles: 1.0},
            tanimoto_cutoff,
        )
        candidates_by_round[0] = seed_table
        best_score = seed_scores[seed_smiles][directing.annotator_id]
        if math.isnan(best_score):
            raise ValueError(
                f"seed_smiles could not be scored by directing annotator "
                f"{directing.annotator_id!r}"
            )
        if not seed_table.iloc[0]["passes_constraints"]:
            logger.warning("Round 0: seed does not satisfy all constraints")
        history.append(
            {
                "round": 0,
                "smiles": seed_smiles,
                "score": best_score,
                "is_new_best": True,
            }
        )
        logger.info(f"Round 0: seed score = {best_score:.3f}")

    for round_num in range(1, n_rounds + 1):
        by_model = generator.generate_by_model(best_smiles)
        source_by_smiles: dict[str, list[str]] = {}
        for model_id, smiles_list in by_model.items():
            for smi in smiles_list:
                source_by_smiles.setdefault(smi, []).append(model_id)
        joined_source = {
            smi: "|".join(sources) for smi, sources in source_by_smiles.items()
        }

        scores = _score_candidates(list(joined_source), annotators)
        tanimoto = (
            tanimoto_to_seed(seed_smiles, list(joined_source))
            if seed_smiles is not None
            else None
        )
        round_table = _round_table(
            joined_source, scores, annotators, tanimoto, tanimoto_cutoff
        )
        candidates_by_round[round_num] = round_table

        if round_table.empty:
            logger.warning(f"Round {round_num}: no candidates generated, stopping.")
            break

        eligible = round_table[
            round_table["passes_constraints"]
            & round_table[directing.annotator_id].notna()
        ]
        if eligible.empty:
            logger.warning(
                f"Round {round_num}: no candidate satisfied all constraints, stopping."
            )
            break

        winner_idx = eligible[directing.annotator_id].idxmax()
        round_best_smiles = eligible.loc[winner_idx, "smiles"]
        round_best_score = eligible.loc[winner_idx, directing.annotator_id]
        is_new_best = round_best_score > best_score + tolerance

        history.append(
            {
                "round": round_num,
                "smiles": round_best_smiles,
                "score": round_best_score,
                "is_new_best": is_new_best,
            }
        )

        if is_new_best:
            logger.success(f"Round {round_num}: improved to {round_best_score:.3f}")
            best_smiles, best_score = round_best_smiles, round_best_score
        else:
            logger.info(f"Round {round_num}: no improvement, stopping.")
            break

    return pd.DataFrame(history), candidates_by_round


def _validate_annotators(annotators: list[AnnotatorSpec]) -> AnnotatorSpec:
    """Check the list-level invariants hill_climb depends on, return the directing spec."""
    ids = [spec.annotator_id for spec in annotators]
    if len(set(ids)) != len(ids):
        raise ValueError("annotators must have unique annotator_id values")

    reserved = _RESERVED_COLUMNS & set(ids)
    if reserved:
        raise ValueError(
            f"annotator_id cannot be one of the reserved names: {reserved}"
        )

    directing = [spec for spec in annotators if spec.role == "directing"]
    if len(directing) != 1:
        raise ValueError(
            f"annotators must have exactly one directing spec, found {len(directing)}"
        )
    return directing[0]


def _score_candidates(
    smiles_list: list[str], annotators: list[AnnotatorSpec]
) -> dict[str, dict[str, float]]:
    """Score every candidate with every annotator, NaN where an annotator drops it."""
    scores_by_id = {
        spec.annotator_id: spec.annotator.score(smiles_list) for spec in annotators
    }
    return {
        smi: {
            aid: scores.get(smi, float("nan")) for aid, scores in scores_by_id.items()
        }
        for smi in smiles_list
    }


def _passes_constraints(
    smi: str,
    scores: dict[str, dict[str, float]],
    controlling: list[AnnotatorSpec],
    tanimoto: dict[str, float] | None,
    tanimoto_cutoff: float | None,
) -> bool:
    """A missing controlling score or missing/low similarity counts as failing."""
    for spec in controlling:
        value = scores[smi].get(spec.annotator_id)
        if value is None or math.isnan(value):
            return False
        if spec.min is not None and value < spec.min:
            return False
        if spec.max is not None and value > spec.max:
            return False
    if tanimoto_cutoff is not None:
        value = None if tanimoto is None else tanimoto.get(smi)
        if value is None or value < tanimoto_cutoff:
            return False
    return True


def _round_table(
    source_by_smiles: dict[str, str],
    scores: dict[str, dict[str, float]],
    annotators: list[AnnotatorSpec],
    tanimoto: dict[str, float] | None,
    tanimoto_cutoff: float | None,
) -> pd.DataFrame:
    """Assemble one round's candidate table: smiles, source, annotator values, verdicts."""
    controlling = [spec for spec in annotators if spec.role == "controlling"]
    rows = []
    for smi, source in source_by_smiles.items():
        row = {"smiles": smi, "source": source, **scores[smi]}
        row["passes_constraints"] = _passes_constraints(
            smi, scores, controlling, tanimoto, tanimoto_cutoff
        )
        if tanimoto is not None:
            row["tanimoto_to_seed"] = tanimoto.get(smi, float("nan"))
        rows.append(row)
    return pd.DataFrame(rows)
