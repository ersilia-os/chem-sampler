import math
import os
from typing import Literal

import pandas as pd

from .models.spec import AnnotatorSpec, Direction
from .utils.logging import logger
from .utils.similarity import tanimoto_to_seed

_RESERVED_COLUMNS = {"smiles", "source", "cutoffs_satisfied", "tanimoto_to_seed"}


def hill_climb(
    generator,
    annotators: list[AnnotatorSpec],
    *,
    mode: Literal["joint", "sequential"],
    seed_smiles: str | None = None,
    n_rounds: int = 5,
    tolerance: float = 0.0,
    tanimoto_cutoff: float | None = None,
) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    """
    Iteratively improve a seed molecule by alternating generation and scoring.

    Every annotator in `annotators` is uniform: a cutoff and a direction, nothing
    more. How they combine into a round's winner is decided by `mode`:

    - "sequential": annotators are optimized one at a time, in list order. Each
      finished stage's achieved value becomes a hard floor for every later stage
      (the achieved value, not just its original cutoff) - a later stage can
      never trade away an earlier gain. A stage whose own cutoff nothing can
      satisfy, not even the entering candidate, stops the whole run.
    - "joint": all annotators optimized together in one search. A candidate's
      score is its count of cutoffs satisfied; ties break on whichever eligible
      candidate is encountered first. `tanimoto_cutoff` is a hard gate in both
      modes but is never counted toward this total.

    Parameters
    ----------
    generator : object
        Any object exposing `.generate_by_model(seed_smiles: str | None) ->
        dict[str, list[str]]` (e.g. `HubGenerator`, `ChemblSampler`, or a
        `GeneratorPool`).
    annotators : list[AnnotatorSpec]
        At least one. In "sequential" mode, list order is the priority order.
    mode : {"joint", "sequential"}
        How annotators combine. No default - with a single annotator the two
        modes diverge sharply (sequential maximizes/minimizes it directly;
        joint only checks pass/fail against its cutoff and picks an eligible
        candidate arbitrarily), so the choice is never silently assumed.
    seed_smiles : str, optional
        SMILES of the starting molecule. If `None`, the first round's winner is
        unconditionally the new best.
    n_rounds : int, optional
        Maximum rounds to run, by default 5. In "sequential" mode this budget
        applies separately to each stage.
    tolerance : float, optional
        Minimum margin a round must beat the current best by to count as an
        improvement, by default 0.0.
    tanimoto_cutoff : float, optional
        Minimum Tanimoto similarity to `seed_smiles` for a candidate to be
        eligible at all. Requires `seed_smiles`.

    Returns
    -------
    summary : pandas.DataFrame
        One row per round attempted, with columns `round`, `smiles`, `score`,
        `is_new_best`, `active_annotator_id` (the annotator that round was
        driving; `None` for every "joint"-mode row).
    candidates_by_round : dict[int, pandas.DataFrame]
        Every candidate considered each round (round 0 is the seed, if given),
        continuously numbered across stage boundaries in "sequential" mode.
        Columns: `smiles`, `source`, one column per `annotators` entry's
        `annotator_id`, `cutoffs_satisfied`, and `tanimoto_to_seed` (only if a
        seed was given; always against the original seed, never the rolling
        best-so-far molecule).

    Raises
    ------
    ValueError
        If `annotators` is empty, has duplicate or reserved `annotator_id`
        values, `mode` isn't "joint"/"sequential", `tanimoto_cutoff` is given
        without a seed, or an entering candidate can't be scored by the
        annotator whose stage it's entering.
    """
    _validate_annotators(annotators)
    if mode not in ("joint", "sequential"):
        raise ValueError(f"mode must be 'joint' or 'sequential', got {mode!r}")
    if tanimoto_cutoff is not None and seed_smiles is None:
        raise ValueError("tanimoto_cutoff requires seed_smiles")

    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    best_row: dict[str, float] = {}
    best_score = float("-inf")

    if seed_smiles is not None:
        seed_scores = _score_candidates([seed_smiles], annotators)
        round0_table = _round_table(
            {seed_smiles: "seed"}, seed_scores, annotators, {seed_smiles: 1.0}
        )
        candidates_by_round[0] = round0_table
        seed_row = round0_table.iloc[0]
        best_row = {
            spec.annotator_id: seed_row[spec.annotator_id] for spec in annotators
        }
        active_id = None if mode == "joint" else annotators[0].annotator_id
        best_score = (
            seed_row["cutoffs_satisfied"] if mode == "joint" else seed_row[active_id]
        )
        history.append(
            {
                "round": 0,
                "smiles": seed_smiles,
                "score": best_score,
                "is_new_best": True,
                "active_annotator_id": active_id,
            }
        )
        logger.info(f"Round 0: seed score = {best_score}")

    if mode == "joint":
        stage_history, stage_candidates = _run_joint(
            generator,
            annotators,
            seed_smiles,
            n_rounds,
            tolerance,
            tanimoto_cutoff,
            best_smiles=seed_smiles,
            best_score=best_score,
            best_row=best_row,
        )
    else:
        stage_history, stage_candidates = _run_sequential(
            generator,
            annotators,
            seed_smiles,
            n_rounds,
            tolerance,
            tanimoto_cutoff,
            best_smiles=seed_smiles,
            best_row=best_row,
        )

    history.extend(stage_history)
    candidates_by_round.update(stage_candidates)
    return pd.DataFrame(history), candidates_by_round


def write_results(
    summary: pd.DataFrame,
    candidates_by_round: dict[int, pd.DataFrame],
    output_dir: str,
) -> None:
    """
    Write a `hill_climb()` result as `summary.csv` plus `round<n>.csv` per round.

    Parameters
    ----------
    summary : pandas.DataFrame
        First element of `hill_climb()`'s return value.
    candidates_by_round : dict[int, pandas.DataFrame]
        Second element of `hill_climb()`'s return value.
    output_dir : str
        Created (including parents) if missing. Files of the same name are
        overwritten; any other files already present are left alone.
    """
    os.makedirs(output_dir, exist_ok=True)
    summary.to_csv(os.path.join(output_dir, "summary.csv"), index=False)
    for round_num, round_table in candidates_by_round.items():
        round_table.to_csv(
            os.path.join(output_dir, f"round{round_num}.csv"), index=False
        )


def _validate_annotators(annotators: list[AnnotatorSpec]) -> None:
    """Check the list-level invariants hill_climb depends on."""
    if not annotators:
        raise ValueError("annotators must not be empty")

    ids = [spec.annotator_id for spec in annotators]
    if len(set(ids)) != len(ids):
        raise ValueError("annotators must have unique annotator_id values")

    reserved = _RESERVED_COLUMNS & set(ids)
    if reserved:
        raise ValueError(
            f"annotator_id cannot be one of the reserved names: {reserved}"
        )


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


def _cutoff_satisfied(value: float, spec: AnnotatorSpec) -> bool:
    """A missing (NaN) value always fails, regardless of how permissive the cutoff is."""
    if value is None or math.isnan(value):
        return False
    if spec.direction == "higher":
        return value >= spec.cutoff
    return value <= spec.cutoff


def _cutoff_satisfied_mask(
    values: pd.Series, cutoff: float, direction: Direction
) -> pd.Series:
    """Vectorized `_cutoff_satisfied`; NaN comparisons are already False in pandas."""
    if direction == "higher":
        return values >= cutoff
    return values <= cutoff


def _cutoffs_satisfied_count(
    scores_row: dict[str, float], annotators: list[AnnotatorSpec]
) -> int:
    return sum(
        _cutoff_satisfied(scores_row[spec.annotator_id], spec) for spec in annotators
    )


def _round_table(
    source_by_smiles: dict[str, str],
    scores: dict[str, dict[str, float]],
    annotators: list[AnnotatorSpec],
    tanimoto: dict[str, float] | None,
) -> pd.DataFrame:
    """Assemble one round's candidate table: smiles, source, annotator values, cutoff count."""
    rows = []
    for smi, source in source_by_smiles.items():
        row = {"smiles": smi, "source": source, **scores[smi]}
        row["cutoffs_satisfied"] = _cutoffs_satisfied_count(scores[smi], annotators)
        if tanimoto is not None:
            row["tanimoto_to_seed"] = tanimoto.get(smi, float("nan"))
        rows.append(row)
    return pd.DataFrame(rows)


def _tanimoto_eligible(
    round_table: pd.DataFrame, tanimoto_cutoff: float | None
) -> pd.Series:
    """All True if no cutoff is given; else Tanimoto-to-seed must clear it (missing -> False)."""
    if tanimoto_cutoff is None:
        return pd.Series(True, index=round_table.index)
    return round_table["tanimoto_to_seed"] >= tanimoto_cutoff


def _improved(
    candidate_score: float,
    best_score: float,
    direction: Direction | None,
    tolerance: float,
) -> bool:
    """direction=None (joint mode's count) and "higher" both maximize; "lower" minimizes."""
    if direction == "lower":
        return candidate_score < best_score - tolerance
    return candidate_score > best_score + tolerance


def _run_rounds(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    start_round: int,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
    pick_winner,
    winner_direction: Direction | None,
    active_annotator_id: str | None,
) -> tuple[list[dict], dict[int, pd.DataFrame], str | None, float, dict[str, float]]:
    """Shared loop: generate -> score -> pick_winner -> compare. Used by both modes."""
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}

    for offset in range(n_rounds):
        round_num = start_round + offset
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
        round_table = _round_table(joined_source, scores, annotators, tanimoto)
        candidates_by_round[round_num] = round_table

        if round_table.empty:
            logger.warning(f"Round {round_num}: no candidates generated, stopping.")
            break

        winner = pick_winner(round_table)
        if winner is None:
            logger.warning(f"Round {round_num}: no eligible candidate, stopping.")
            break
        round_best_smiles, round_best_score, round_best_row = winner

        is_new_best = _improved(
            round_best_score, best_score, winner_direction, tolerance
        )
        history.append(
            {
                "round": round_num,
                "smiles": round_best_smiles,
                "score": round_best_score,
                "is_new_best": is_new_best,
                "active_annotator_id": active_annotator_id,
            }
        )

        if is_new_best:
            logger.success(f"Round {round_num}: improved to {round_best_score}")
            best_smiles, best_score, best_row = (
                round_best_smiles,
                round_best_score,
                round_best_row,
            )
        else:
            logger.info(f"Round {round_num}: no improvement, stopping.")
            break

    return history, candidates_by_round, best_smiles, best_score, best_row


def _row_values(row: pd.Series, annotators: list[AnnotatorSpec]) -> dict[str, float]:
    """Every annotator's value for one round-table row, as a plain dict."""
    return {spec.annotator_id: row[spec.annotator_id] for spec in annotators}


def _run_joint(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    def pick_winner(round_table: pd.DataFrame):
        eligible = round_table[_tanimoto_eligible(round_table, tanimoto_cutoff)]
        if eligible.empty:
            return None
        idx = eligible["cutoffs_satisfied"].idxmax()
        row = eligible.loc[idx]
        return row["smiles"], row["cutoffs_satisfied"], _row_values(row, annotators)

    history, candidates_by_round, _, _, _ = _run_rounds(
        generator,
        annotators,
        seed_smiles,
        n_rounds,
        tolerance,
        start_round=1,
        best_smiles=best_smiles,
        best_score=best_score,
        best_row=best_row,
        pick_winner=pick_winner,
        winner_direction=None,
        active_annotator_id=None,
    )
    return history, candidates_by_round


def _run_sequential(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    best_smiles: str | None,
    best_row: dict[str, float],
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    floors: list[tuple[str, float, Direction]] = []
    next_round = 1

    for spec in annotators:
        if best_smiles is None:
            best_score = float("-inf") if spec.direction == "higher" else float("inf")
            entry_ok = True  # no entry candidate to check; the "no candidates" path handles failure
        else:
            entry_value = best_row.get(spec.annotator_id, float("nan"))
            if math.isnan(entry_value):
                raise ValueError(
                    f"{best_smiles!r} could not be scored by {spec.annotator_id!r} "
                    "entering its stage"
                )
            best_score = entry_value
            entry_ok = _cutoff_satisfied(entry_value, spec)

        floors_snapshot = list(floors)

        def pick_winner(
            round_table: pd.DataFrame, spec=spec, floors_snapshot=floors_snapshot
        ):
            mask = _tanimoto_eligible(round_table, tanimoto_cutoff)
            mask &= _cutoff_satisfied_mask(
                round_table[spec.annotator_id], spec.cutoff, spec.direction
            )
            for floor_id, floor_value, floor_direction in floors_snapshot:
                mask &= _cutoff_satisfied_mask(
                    round_table[floor_id], floor_value, floor_direction
                )
            eligible = round_table[mask]
            if eligible.empty:
                return None
            idx = (
                eligible[spec.annotator_id].idxmax()
                if spec.direction == "higher"
                else eligible[spec.annotator_id].idxmin()
            )
            row = eligible.loc[idx]
            return row["smiles"], row[spec.annotator_id], _row_values(row, annotators)

        (
            stage_history,
            stage_candidates,
            new_best_smiles,
            new_best_score,
            new_best_row,
        ) = _run_rounds(
            generator,
            annotators,
            seed_smiles,
            n_rounds,
            tolerance,
            start_round=next_round,
            best_smiles=best_smiles,
            best_score=best_score,
            best_row=best_row,
            pick_winner=pick_winner,
            winner_direction=spec.direction,
            active_annotator_id=spec.annotator_id,
        )
        history.extend(stage_history)
        candidates_by_round.update(stage_candidates)
        next_round += len(stage_candidates)

        if new_best_smiles is None:
            logger.warning(
                f"{spec.annotator_id}: no candidate found, stopping the whole run."
            )
            break

        if new_best_smiles == best_smiles:
            if not entry_ok:
                logger.warning(
                    f"{spec.annotator_id}: no candidate satisfies its cutoff, "
                    "stopping the whole run."
                )
                break
            achieved = best_score
        else:
            achieved = new_best_score
            best_smiles, best_row = new_best_smiles, new_best_row

        floors.append((spec.annotator_id, achieved, spec.direction))

    return history, candidates_by_round
