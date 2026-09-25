import logging
import math
import os
from typing import Literal

import pandas as pd

from .models.spec import AnnotatorSpec, Direction
from .utils.logging import logger
from .utils.similarity import tanimoto_to_seed

_RESERVED_COLUMNS = {
    "smiles",
    "source",
    "cutoffs_satisfied",
    "weighted_score",
    "tanimoto_to_seed",
    "tanimoto_to_original_seed",
}


def hill_climb(
    generator,
    annotators: list[AnnotatorSpec],
    *,
    mode: Literal["joint", "sequential", "weighted"],
    seed_smiles: str | None = None,
    original_seed_smiles: str | None = None,
    n_rounds: int = 5,
    tolerance: float = 0.0,
    tanimoto_cutoff: float | None = None,
    tanimoto_direction: Direction = "higher",
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
      candidate is encountered first.
    - "weighted": all annotators optimized together as one scalarized score,
      `sum(spec.weight * value)` per annotator, with `value` negated first for
      any `direction="lower"` annotator - so maximizing this sum always pushes
      every annotator toward its own better direction, regardless of mix.
      `weight` comes from `AnnotatorSpec.weight` (default 1.0), and is otherwise
      unused by every mode. No cross-annotator normalization is applied, so
      weight magnitudes should be chosen relative to each annotator's own scale
      (e.g. a 0-1 QED score needs a much larger weight than a score in the
      hundreds to have comparable pull). A candidate any weighted annotator
      can't score is ineligible to win that round (its total is NaN).

    `tanimoto_cutoff` is a hard gate in every mode, in the direction set by
    `tanimoto_direction`, but is never counted toward any mode's own score.

    Parameters
    ----------
    generator : object
        Any object exposing `.generate_by_model(seed_smiles: str | None) ->
        dict[str, list[str]]` (e.g. `HubGenerator`, `ChemblSampler`, or a
        `GeneratorPool`).
    annotators : list[AnnotatorSpec]
        At least one. In "sequential" mode, list order is the priority order.
    mode : {"joint", "sequential", "weighted"}
        How annotators combine. No default - with a single annotator the modes
        can still diverge sharply (sequential/weighted maximize or minimize it
        directly; joint only checks pass/fail against its cutoff and picks an
        eligible candidate arbitrarily), so the choice is never silently assumed.
    seed_smiles : str, optional
        SMILES of the starting molecule. If `None`, the first round's winner is
        unconditionally the new best.
    original_seed_smiles : str, optional
        SMILES of the true original molecule, for tracking similarity across a
        manually re-seeded chain of `hill_climb` calls. Adds
        `tanimoto_to_original_seed` to `candidates_by_round`, independent of
        `seed_smiles`/`tanimoto_to_seed`; never affects eligibility or the
        search itself.
    n_rounds : int, optional
        Maximum rounds to run, by default 5. In "sequential" mode this budget
        applies separately to each stage.
    tolerance : float, optional
        Minimum margin a round must beat the current best by to count as an
        improvement, by default 0.0.
    tanimoto_cutoff : float, optional
        Tanimoto similarity to `seed_smiles` a candidate must clear to be
        eligible at all, gated by `tanimoto_direction`. Requires `seed_smiles`.
    tanimoto_direction : {"higher", "lower"}, optional
        "higher" (default) keeps candidates at least as similar to `seed_smiles`
        as `tanimoto_cutoff` - the pre-existing similarity-floor behavior.
        "lower" keeps candidates at most that similar instead, pushing the
        search toward novelty. Ignored if `tanimoto_cutoff` is `None`.

    Returns
    -------
    summary : pandas.DataFrame
        One row per round attempted, with columns `round`, `smiles`, `score`,
        `is_new_best`, `active_annotator_id` (the annotator that round was
        driving; `None` for every "joint"/"weighted"-mode row).
    candidates_by_round : dict[int, pandas.DataFrame]
        Every candidate considered each round (round 0 is the seed, if given),
        continuously numbered across stage boundaries in "sequential" mode.
        Columns: `smiles`, `source`, one column per `annotators` entry's
        `annotator_id`, `cutoffs_satisfied`, `weighted_score` (both always
        present, regardless of `mode`), `tanimoto_to_seed` (only if a seed was
        given; always against `seed_smiles` itself, never the rolling
        best-so-far molecule), and `tanimoto_to_original_seed` (only if
        `original_seed_smiles` was given).

    Raises
    ------
    ValueError
        If `annotators` is empty, has duplicate or reserved `annotator_id`
        values, `mode` isn't "joint"/"sequential"/"weighted", `tanimoto_direction`
        isn't "higher"/"lower", `tanimoto_cutoff` is given without a seed, or an
        entering candidate can't be scored by the annotator whose stage it's
        entering.
    """
    _validate_annotators(annotators)
    if mode not in ("joint", "sequential", "weighted"):
        raise ValueError(
            f"mode must be 'joint', 'sequential', or 'weighted', got {mode!r}"
        )
    if tanimoto_direction not in ("higher", "lower"):
        raise ValueError(
            f"tanimoto_direction must be 'higher' or 'lower', got {tanimoto_direction!r}"
        )
    if tanimoto_cutoff is not None and seed_smiles is None:
        raise ValueError("tanimoto_cutoff requires seed_smiles")

    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    best_row: dict[str, float] = {}
    best_score = float("-inf")

    if seed_smiles is not None:
        seed_scores = _score_candidates([seed_smiles], annotators)
        original_seed_tanimoto = (
            tanimoto_to_seed(original_seed_smiles, [seed_smiles])
            if original_seed_smiles is not None
            else None
        )
        round0_table = _round_table(
            {seed_smiles: "seed"},
            seed_scores,
            annotators,
            {seed_smiles: 1.0},
            original_seed_tanimoto,
        )
        candidates_by_round[0] = round0_table
        seed_row = round0_table.iloc[0]
        best_row = {
            spec.annotator_id: seed_row[spec.annotator_id] for spec in annotators
        }
        active_id = (
            None if mode in ("joint", "weighted") else annotators[0].annotator_id
        )
        if mode == "joint":
            best_score = seed_row["cutoffs_satisfied"]
        elif mode == "weighted":
            best_score = seed_row["weighted_score"]
        else:
            best_score = seed_row[active_id]
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
            original_seed_smiles,
            n_rounds,
            tolerance,
            tanimoto_cutoff,
            tanimoto_direction,
            best_smiles=seed_smiles,
            best_score=best_score,
            best_row=best_row,
        )
    elif mode == "weighted":
        stage_history, stage_candidates = _run_weighted(
            generator,
            annotators,
            seed_smiles,
            original_seed_smiles,
            n_rounds,
            tolerance,
            tanimoto_cutoff,
            tanimoto_direction,
            best_smiles=seed_smiles,
            best_score=best_score,
            best_row=best_row,
        )
    else:
        stage_history, stage_candidates = _run_sequential(
            generator,
            annotators,
            seed_smiles,
            original_seed_smiles,
            n_rounds,
            tolerance,
            tanimoto_cutoff,
            tanimoto_direction,
            best_smiles=seed_smiles,
            best_row=best_row,
        )

    history.extend(stage_history)
    candidates_by_round.update(stage_candidates)
    return pd.DataFrame(history), candidates_by_round


def annotate(smiles: str, annotators: list[AnnotatorSpec]) -> pd.DataFrame:
    """
    Score one molecule against a set of annotators.

    A fast, single-shot counterpart to `hill_climb`: no generation, no rounds.
    Useful for checking a candidate seed's baseline, or spot-checking any
    molecule, before committing to a generator and a full run.

    Parameters
    ----------
    smiles : str
        SMILES of the molecule to score.
    annotators : list[AnnotatorSpec]
        At least one. See `hill_climb` for the shared invariants (unique,
        non-reserved `annotator_id` values).

    Returns
    -------
    pandas.DataFrame
        Exactly one row, shaped like a `hill_climb` round table: `smiles`,
        `source` (always "input"), one column per `annotators` entry's
        `annotator_id` (`NaN` where that annotator can't score `smiles`), and
        `cutoffs_satisfied`.

    Raises
    ------
    ValueError
        If `annotators` is empty or has duplicate/reserved `annotator_id`
        values.
    """
    _validate_annotators(annotators)
    scores = _score_candidates([smiles], annotators)
    return _round_table({smiles: "input"}, scores, annotators, tanimoto=None)


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


def write_results_incremental(
    summary: list[dict],
    candidates_by_round: dict[int, pd.DataFrame],
    output_dir: str,
) -> None:
    """
    Write partial `hill_climb()` results incrementally during a run.

    Use after each stage completes to save progress; call `write_results()` at the
    end with the final DataFrame. This allows partial results to survive interruption.

    Parameters
    ----------
    summary : list[dict]
        In-progress history list from `hill_climb()`.
    candidates_by_round : dict[int, pandas.DataFrame]
        Candidates collected so far.
    output_dir : str
        Same as `write_results()`.
    """
    if not summary:
        return
    os.makedirs(output_dir, exist_ok=True)
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(output_dir, "summary.csv"), index=False)
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


def _weighted_score(
    scores_row: dict[str, float], annotators: list[AnnotatorSpec]
) -> float:
    """Sum of weight * signed value ("lower"-direction values are negated first).

    A NaN from any annotator propagates to the whole sum, which sorts and
    compares last everywhere else in this module - the same "missing always
    fails" rule `_cutoff_satisfied` already applies, with no special-casing
    needed here.
    """
    total = 0.0
    for spec in annotators:
        value = scores_row[spec.annotator_id]
        signed = value if spec.direction == "higher" else -value
        total += spec.weight * signed
    return total


def _round_table(
    source_by_smiles: dict[str, str],
    scores: dict[str, dict[str, float]],
    annotators: list[AnnotatorSpec],
    tanimoto: dict[str, float] | None,
    original_seed_tanimoto: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Assemble one round's candidate table: smiles, source, annotator values, cutoff count."""
    rows = []
    for smi, source in source_by_smiles.items():
        row = {"smiles": smi, "source": source, **scores[smi]}
        row["cutoffs_satisfied"] = _cutoffs_satisfied_count(scores[smi], annotators)
        row["weighted_score"] = _weighted_score(scores[smi], annotators)
        if tanimoto is not None:
            row["tanimoto_to_seed"] = tanimoto.get(smi, float("nan"))
        if original_seed_tanimoto is not None:
            row["tanimoto_to_original_seed"] = original_seed_tanimoto.get(
                smi, float("nan")
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _sort_round_table(
    round_table: pd.DataFrame,
    active_annotator_id: str | None,
    winner_direction: Direction | None,
    sort_column: str = "cutoffs_satisfied",
) -> pd.DataFrame:
    """Order a round's candidates best-first, by whatever that round is optimizing.

    Sequential mode (`active_annotator_id` given): that annotator's own column,
    ascending if its direction is "lower", else descending. Joint/weighted mode
    (`active_annotator_id` is None): `sort_column` descending - `cutoffs_satisfied`
    for joint mode, `weighted_score` for weighted mode; the same column that
    mode's own `pick_winner` already uses via `idxmax`. NaN in the sort column
    always sorts last (pandas default), so a candidate the active annotator
    couldn't score is never treated as "best".
    """
    if round_table.empty:
        return round_table
    if active_annotator_id is not None:
        column, ascending = active_annotator_id, winner_direction == "lower"
    else:
        column, ascending = sort_column, False
    return round_table.sort_values(
        column, ascending=ascending, kind="stable"
    ).reset_index(drop=True)


def _tanimoto_eligible(
    round_table: pd.DataFrame,
    tanimoto_cutoff: float | None,
    tanimoto_direction: Direction,
) -> pd.Series:
    """All True if no cutoff is given; else Tanimoto-to-seed must clear it (missing -> False)."""
    if tanimoto_cutoff is None:
        return pd.Series(True, index=round_table.index)
    return _cutoff_satisfied_mask(
        round_table["tanimoto_to_seed"], tanimoto_cutoff, tanimoto_direction
    )


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
    original_seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    start_round: int,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
    pick_winner,
    winner_direction: Direction | None,
    active_annotator_id: str | None,
    active_cutoff: float | None = None,
    sort_column: str = "cutoffs_satisfied",
) -> tuple[list[dict], dict[int, pd.DataFrame], str | None, float, dict[str, float]]:
    """Shared loop: generate -> score -> pick_winner -> compare. Used by every mode."""
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}

    for offset in range(n_rounds):
        round_num = start_round + offset

        if logger._quiet_mode:
            logger.logger.setLevel(logging.ERROR)
        by_model = generator.generate_by_model(best_smiles)
        if logger._quiet_mode:
            logger.logger.setLevel(logging.INFO)
            logger.info(f"Round {round_num}: running...")
        else:
            if hasattr(generator, "generators"):
                gen_ids = [g.model_id for g in generator.generators]
            else:
                gen_ids = list(by_model.keys()) if by_model else ["unknown"]
            gen_str = ", ".join(gen_ids)
            ann_ids = [spec.annotator_id for spec in annotators]
            ann_str = ", ".join(ann_ids)
            logger.info(
                f"Round {round_num}: Generators: {gen_str}. Annotators: {ann_str}. Running..."
            )

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
        original_seed_tanimoto = (
            tanimoto_to_seed(original_seed_smiles, list(joined_source))
            if original_seed_smiles is not None
            else None
        )
        round_table = _round_table(
            joined_source, scores, annotators, tanimoto, original_seed_tanimoto
        )
        round_table = _sort_round_table(
            round_table, active_annotator_id, winner_direction, sort_column
        )
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
            n_candidates = len(round_table)
            if logger._quiet_mode:
                logger.success(
                    f"Round {round_num}: improved to {round_best_score:.2f} ({n_candidates} candidates)"
                )
            else:
                logger.success(f"Round {round_num}: improved to {round_best_score}")
            best_smiles, best_score, best_row = (
                round_best_smiles,
                round_best_score,
                round_best_row,
            )

            if (
                active_cutoff is not None
                and winner_direction is not None
                and math.isfinite(active_cutoff)
            ):
                cutoff_met = (
                    round_best_score >= active_cutoff
                    if winner_direction == "higher"
                    else round_best_score <= active_cutoff
                )
                if cutoff_met:
                    logger.info(
                        f"{active_annotator_id}: cutoff {active_cutoff} satisfied, "
                        "stopping this stage."
                    )
                    break
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
    original_seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: Direction,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    def pick_winner(round_table: pd.DataFrame):
        eligible = round_table[
            _tanimoto_eligible(round_table, tanimoto_cutoff, tanimoto_direction)
        ]
        if eligible.empty:
            return None
        idx = eligible["cutoffs_satisfied"].idxmax()
        row = eligible.loc[idx]
        return row["smiles"], row["cutoffs_satisfied"], _row_values(row, annotators)

    history, candidates_by_round, _, _, _ = _run_rounds(
        generator,
        annotators,
        seed_smiles,
        original_seed_smiles,
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


def _run_weighted(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: Direction,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    def pick_winner(round_table: pd.DataFrame):
        mask = _tanimoto_eligible(round_table, tanimoto_cutoff, tanimoto_direction)
        mask &= round_table["weighted_score"].notna()
        eligible = round_table[mask]
        if eligible.empty:
            return None
        idx = eligible["weighted_score"].idxmax()
        row = eligible.loc[idx]
        return row["smiles"], row["weighted_score"], _row_values(row, annotators)

    history, candidates_by_round, _, _, _ = _run_rounds(
        generator,
        annotators,
        seed_smiles,
        original_seed_smiles,
        n_rounds,
        tolerance,
        start_round=1,
        best_smiles=best_smiles,
        best_score=best_score,
        best_row=best_row,
        pick_winner=pick_winner,
        winner_direction=None,
        active_annotator_id=None,
        sort_column="weighted_score",
    )
    return history, candidates_by_round


def _run_sequential(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: Direction,
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
            mask = _tanimoto_eligible(round_table, tanimoto_cutoff, tanimoto_direction)
            mask &= _cutoff_satisfied_mask(
                round_table[spec.annotator_id], spec.cutoff, spec.direction
            )
            for floor_id, floor_value, floor_direction in floors_snapshot:
                mask &= _cutoff_satisfied_mask(
                    round_table[floor_id], floor_value, floor_direction
                )
            eligible = round_table[mask]
            if eligible.empty:
                # No candidate clears this stage's own cutoff even with floors
                # intact. Still report an improvement if one exists that respects
                # every prior floor - floors are a hard invariant and are never
                # dropped, even in this fallback (own cutoff is still excluded;
                # that's the whole point of the fallback).
                floor_mask = _tanimoto_eligible(
                    round_table, tanimoto_cutoff, tanimoto_direction
                )
                for floor_id, floor_value, floor_direction in floors_snapshot:
                    floor_mask &= _cutoff_satisfied_mask(
                        round_table[floor_id], floor_value, floor_direction
                    )
                if floor_mask.any():
                    best_overall = round_table[floor_mask]
                    idx = (
                        best_overall[spec.annotator_id].idxmax()
                        if spec.direction == "higher"
                        else best_overall[spec.annotator_id].idxmin()
                    )
                    row = best_overall.loc[idx]
                    if logger._quiet_mode:
                        logger.warning(
                            f"{spec.annotator_id}: {row[spec.annotator_id]:.2f} (cutoff {spec.cutoff} not met)"
                        )
                    else:
                        logger.warning(
                            f"{spec.annotator_id}: improved to {row[spec.annotator_id]}, "
                            f"but did not satisfy cutoff of {spec.cutoff} (direction: {spec.direction}). "
                            "Continuing to next round."
                        )
                    return (
                        row["smiles"],
                        row[spec.annotator_id],
                        _row_values(row, annotators),
                    )
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
            original_seed_smiles,
            n_rounds,
            tolerance,
            start_round=next_round,
            best_smiles=best_smiles,
            best_score=best_score,
            best_row=best_row,
            pick_winner=pick_winner,
            winner_direction=spec.direction,
            active_annotator_id=spec.annotator_id,
            active_cutoff=spec.cutoff,
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
