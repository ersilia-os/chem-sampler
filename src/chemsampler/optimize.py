import dataclasses
import logging
import math
import os
from functools import partial

import pandas as pd

from . import execution
from .models.annotator import NEEDS_SEED, TanimotoAnnotator
from .models.generator import SeedRequiredError
from .models.similarity import HubSimilarity
from .models.spec import (
    AnnotatorDirection,
    AnnotatorSpec,
    Cutoff,
    Direction,
    has_cutoff,
)
from .utils.logging import logger
from .utils.similarity import tanimoto_to_seed

_RESERVED_COLUMNS = {
    "smiles",
    "source",
    "cutoff_satisfied",
    "guardrails_satisfied",
    "tanimoto_to_seed",
}

#: What stops a run midway: Ctrl-C, or a model that fails. The rounds finished
#: before it travel on the exception, so the caller can still save them.
_STOPS = (KeyboardInterrupt, RuntimeError, ValueError)


def _keep_finished(error: BaseException, history: list, candidates_by_round: dict):
    """Put the rounds a caller finished before those its callee put on `error`."""
    later_history, later_candidates = getattr(error, "_finished_rounds", ([], {}))
    error._finished_rounds = (
        [*history, *later_history],
        {**candidates_by_round, **later_candidates},
    )


def _table_columns(annotators: list[AnnotatorSpec], tanimoto: bool) -> set[str]:
    """The columns `_round_table` gives a non-empty round."""
    columns = {"smiles", "source", "cutoff_satisfied", "guardrails_satisfied"}
    columns |= {spec.annotator_id for spec in annotators}
    if tanimoto:
        columns.add("tanimoto_to_seed")
    return columns


def _saved_round(checkpoint, round_num: int, columns: set[str]) -> pd.DataFrame | None:
    """
    The table a stopped run saved for this round, or None when it has none.

    Read back exactly as it was made (every float to the last digit, only empty
    cells as missing), so the round's winner and everything after it come out
    as they did the first time.
    """
    if checkpoint is None:
        return None
    path = os.path.join(checkpoint, f"round{round_num}.csv")
    if not os.path.isfile(path):
        return None
    try:
        table = pd.read_csv(
            path,
            float_precision="round_trip",
            keep_default_na=False,
            na_values=[""],
            dtype={"smiles": str, "source": str},
        )
    except pd.errors.EmptyDataError:  # a round that found no candidates
        return pd.DataFrame()
    if set(table.columns) != columns:
        raise ValueError(
            f"{path} was not made by this run's annotators and options: it has "
            f"columns {sorted(table.columns)}, this run makes {sorted(columns)}"
        )
    return table


def _save_round(checkpoint, round_num: int, table: pd.DataFrame) -> None:
    """Write a round's table whole, then rename it: a killed run leaves no half file."""
    if checkpoint is None:
        return
    os.makedirs(checkpoint, exist_ok=True)
    path = os.path.join(checkpoint, f"round{round_num}.csv")
    table.to_csv(path + ".partial", index=False)
    os.replace(path + ".partial", path)


def hill_climb(
    generator,
    annotators: list[AnnotatorSpec],
    *,
    similarity: HubSimilarity | list[HubSimilarity] | None = None,
    seed_smiles: str | None = None,
    n_rounds: int = 5,
    tolerance: float = 0.0,
    checkpoint: str | os.PathLike | None = None,
) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    """
    Iteratively improve a seed molecule by alternating generation and scoring.

    Exactly one annotator is *primary*: the one optimized. Every other one is a
    *guard-rail*, a hard constraint: a candidate that does not satisfy its cutoff
    or range is never eligible to win, in any round, and it takes no part in the
    optimization itself.

    Each round, the candidates proposed around the current best molecule are
    scored; among those that satisfy every guard-rail, the best primary value
    (highest for direction "higher", lowest for "lower") wins. It becomes the
    new best, and the next round's starting molecule, if it beats the current
    best by more than `tolerance`. The run stops when a round does not improve,
    when no candidate is eligible, when the primary's cutoff is met, or after
    `n_rounds`. A primary whose cutoff is NaN has no cutoff: it is optimized
    until a round does not improve or the rounds run out. An infinite cutoff
    never stops a run early either.

    A `TanimotoAnnotator` (`rdkit,tanimoto`) scores each candidate's similarity
    to `seed_smiles` and works like any other annotator: as a guard-rail it
    keeps candidates near the seed (or away from it), as the primary it is
    optimized.

    Parameters
    ----------
    generator : object or None
        Any object exposing `.generate_by_model(seed_smiles: str | None) ->
        dict[str, list[str]]` (e.g. `HubGenerator` or a `GeneratorPool`). May be
        `None` when `similarity` is given.
    annotators : list[AnnotatorSpec]
        Exactly one primary; the rest are guard-rails.
    similarity : HubSimilarity or list of HubSimilarity, optional
        Similarity-search models whose nearest neighbours of the current best
        molecule join each round's candidates next to the generator's, or are
        the only candidates when `generator` is `None`. Kept separate from
        `generator`, which only takes generative models. Skipped with a warning
        when there is no seed to search from.
    seed_smiles : str, optional
        SMILES of the starting molecule. If `None`, the first round's winner is
        unconditionally the new best. Required if a `TanimotoAnnotator` is among
        `annotators`.
    n_rounds : int, optional
        Maximum rounds to run, by default 5.
    tolerance : float, optional
        Minimum margin a round must beat the current best by to count as an
        improvement, by default 0.0.
    checkpoint : str or path-like, optional
        A folder where each round's table is written as `round<n>.csv` as soon
        as it is made, and from which rounds already there are read back instead
        of being run again. A run stopped midway, even killed, continues where
        it stopped when called again with the same arguments and folder (see
        Notes).

    Returns
    -------
    summary : pandas.DataFrame
        One row per round attempted, with columns `round`, `smiles` (the round's
        winner), `score` (its primary value) and `is_new_best`.
    candidates_by_round : dict[int, pandas.DataFrame]
        Every candidate considered each round (round 0 is the seed, if given),
        best first. Columns: `smiles`, `source`, one column per `annotators`
        entry's `annotator_id`, `cutoff_satisfied` (the primary's; empty when it
        has no cutoff), `guardrails_satisfied` (True when every guard-rail is
        satisfied; only such candidates can win) and `tanimoto_to_seed` (only if
        a seed was given; always against `seed_smiles` itself, never the rolling
        best-so-far molecule).

    Raises
    ------
    ValueError
        If neither `generator` nor `similarity` is given, `annotators` is empty,
        has no primary or more than one, has duplicate or reserved
        `annotator_id` values, a `TanimotoAnnotator` is given without a seed, or
        the seed can't be scored by the primary.
    TypeError
        If `generator` is a `HubSimilarity`: pass it as `similarity` instead.

    Notes
    -----
    A run stopped midway by Ctrl-C (`KeyboardInterrupt`) or a failing model
    (`RuntimeError`, `ValueError`) re-raises that exception with the rounds
    finished before it as `partial`: the `(summary, candidates_by_round)` pair
    this function would have returned. It is absent when no round finished.

    With `checkpoint`, nothing but the round tables needs saving to continue: a
    round's winner and when to stop are worked out from the tables, so rounds
    read back lead to exactly the state the first call reached. The round that
    was interrupted is run again whole. Calling again with other arguments than
    the first time gives a run that neither call describes.
    """
    if isinstance(generator, HubSimilarity):
        raise TypeError(
            f"{generator.model_id} is a similarity-search model, not a generator; "
            "pass it as hill_climb(similarity=...)"
        )
    if similarity is None:
        similarity = []
    elif not isinstance(similarity, (list, tuple)):
        similarity = [similarity]
    if generator is None and not similarity:
        raise ValueError("give a generator, a similarity model, or both")
    _validate_annotators(annotators)
    annotators = _with_seed(annotators, seed_smiles)
    primary = _primary(annotators)

    if similarity:
        generator = _CandidateSources(generator, similarity)

    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    best_row: dict[str, float] = {}

    try:
        if seed_smiles is not None:
            columns = _table_columns(annotators, True)
            round0_table = _saved_round(checkpoint, 0, columns)
            if round0_table is not None:
                logger.info("Round 0: the seed's scores, from round0.csv")
            else:
                with execution.stage("Round 0 · annotation"):
                    seed_scores = _score_candidates([seed_smiles], annotators)
                round0_table = _round_table(
                    {seed_smiles: "seed"},
                    seed_scores,
                    annotators,
                    {seed_smiles: 1.0},
                )
                _save_round(checkpoint, 0, round0_table)
            candidates_by_round[0] = round0_table
            seed_row = round0_table.iloc[0]
            failing = [
                spec.annotator_id
                for spec in _guards(annotators)
                if not _cutoff_satisfied(seed_row[spec.annotator_id], spec)
            ]
            if failing:
                logger.warning(
                    f"The seed does not satisfy the guard-rail(s) {', '.join(failing)}; "
                    "its candidates must."
                )
            best_row = {
                spec.annotator_id: seed_row[spec.annotator_id] for spec in annotators
            }
            seed_score = seed_row[primary.annotator_id]
            history.append(
                {
                    "round": 0,
                    "smiles": seed_smiles,
                    "score": seed_score,
                    "is_new_best": True,
                }
            )
            logger.info(f"Round 0: seed score = {seed_score}")

        stage_history, stage_candidates = _run(
            generator,
            annotators,
            seed_smiles,
            n_rounds,
            tolerance,
            best_smiles=seed_smiles,
            best_row=best_row,
            checkpoint=checkpoint,
        )

    except _STOPS as error:
        _keep_finished(error, history, candidates_by_round)
        finished_history, finished_candidates = error._finished_rounds
        if finished_history:
            error.partial = (pd.DataFrame(finished_history), finished_candidates)
        raise

    history.extend(stage_history)
    candidates_by_round.update(stage_candidates)
    return pd.DataFrame(history), candidates_by_round


def annotate(
    smiles: str, annotators: list[AnnotatorSpec], seed_smiles: str | None = None
) -> pd.DataFrame:
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
        At least one, with at most one primary. See `hill_climb` for the shared
        invariants (unique, non-reserved `annotator_id` values).
    seed_smiles : str, optional
        The seed a `TanimotoAnnotator` compares `smiles` with.

    Returns
    -------
    pandas.DataFrame
        Exactly one row, shaped like a `hill_climb` round table: `smiles`,
        `source` (always "input"), one column per `annotators` entry's
        `annotator_id` (`NaN` where that annotator can't score `smiles`),
        `cutoff_satisfied` (the primary's; empty without a primary or a cutoff)
        and `guardrails_satisfied`.

    Raises
    ------
    ValueError
        If `annotators` is empty, has more than one primary or duplicate/reserved
        `annotator_id` values, or has a `TanimotoAnnotator` and no `seed_smiles`
        is given.
    """
    _validate_annotators(annotators, need_primary=False)
    annotators = _with_seed(annotators, seed_smiles)
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


class _CandidateSources:
    """A generator and similarity models behind `hill_climb`'s one-generator interface."""

    def __init__(self, generator, similarity: list[HubSimilarity]):
        self._generator = generator
        self._similarity = similarity
        # The per-round log line reads `.generators`; a lone generator has none.
        generators = (
            [] if generator is None else getattr(generator, "generators", [generator])
        )
        self.generators = [*generators, *similarity]

    def generate_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        # the generators and the similarity search run side by side on a cluster
        calls = [] if self._generator is None else [self._generator.generate_by_model]
        calls.extend(
            partial(self._search_by_model, model) for model in self._similarity
        )
        by_model = {}
        for found in execution.current().map(lambda call: call(seed_smiles), calls):
            by_model.update(found)
        return by_model

    @staticmethod
    def _search_by_model(model, seed_smiles: str | None) -> dict[str, list[str]]:
        try:
            return model.search_by_model(seed_smiles)
        except SeedRequiredError:
            logger.warning(f"{model.model_id}: requires a seed molecule, skipped")
            return {}


def _primaries(annotators: list[AnnotatorSpec]) -> list[AnnotatorSpec]:
    return [spec for spec in annotators if spec.priority == "primary"]


def _primary(annotators: list[AnnotatorSpec]) -> AnnotatorSpec | None:
    """The one annotator a run optimizes, or None (only `annotate` allows none)."""
    primaries = _primaries(annotators)
    return primaries[0] if primaries else None


def _guards(annotators: list[AnnotatorSpec]) -> list[AnnotatorSpec]:
    """The guard-rails: constraints every candidate must satisfy."""
    return [spec for spec in annotators if spec.priority == "guard-rail"]


def _validate_annotators(
    annotators: list[AnnotatorSpec], need_primary: bool = True
) -> None:
    """Check the list-level invariants hill_climb depends on."""
    if not annotators:
        raise ValueError("annotators must not be empty")
    primaries = [spec.annotator_id for spec in _primaries(annotators)]
    if need_primary and not primaries:
        raise ValueError(
            "a primary annotator is required: guard-rails are not optimized"
        )
    if len(primaries) > 1:
        raise ValueError(
            f"one primary annotator only, got {len(primaries)} "
            f"({', '.join(primaries)}): make the others guard-rails"
        )

    ids = [spec.annotator_id for spec in annotators]
    if len(set(ids)) != len(ids):
        raise ValueError("annotators must have unique annotator_id values")

    reserved = _RESERVED_COLUMNS & set(ids)
    if reserved:
        raise ValueError(
            f"annotator_id cannot be one of the reserved names: {reserved}"
        )


def _with_seed(
    annotators: list[AnnotatorSpec], seed_smiles: str | None
) -> list[AnnotatorSpec]:
    """The annotators, with every Tanimoto one comparing candidates to the seed."""
    tanimoto = [s for s in annotators if isinstance(s.annotator, TanimotoAnnotator)]
    if not tanimoto:
        return annotators
    if seed_smiles is None:
        raise ValueError(NEEDS_SEED)
    primary = _primary(annotators)
    if primary in tanimoto and primary.direction == "higher":
        logger.warning(
            f"{primary.annotator_id} (higher) is the primary: nothing can be more "
            "similar to the seed than the seed itself (1.0), so the run cannot "
            "improve on it"
        )
    return [
        dataclasses.replace(s, annotator=s.annotator.against(seed_smiles))
        if s in tanimoto
        else s
        for s in annotators
    ]


def _score_candidates(
    smiles_list: list[str], annotators: list[AnnotatorSpec]
) -> dict[str, dict[str, float]]:
    """Score every candidate with every annotator, NaN where an annotator drops it."""
    # one call per annotator, side by side on a cluster
    scores = execution.current().map(
        lambda spec: spec.annotator.score(smiles_list), annotators
    )
    scores_by_id = {spec.annotator_id: found for spec, found in zip(annotators, scores)}
    return {
        smi: {
            aid: scores.get(smi, float("nan")) for aid, scores in scores_by_id.items()
        }
        for smi in smiles_list
    }


def _value_satisfies(
    value: float, cutoff: Cutoff, direction: AnnotatorDirection
) -> bool:
    """A missing (NaN) value always fails, regardless of how permissive the cutoff is."""
    if value is None or math.isnan(value):
        return False
    if direction == "range":
        low, high = cutoff
        return low <= value <= high
    if direction == "higher":
        return value >= cutoff
    return value <= cutoff


def _cutoff_satisfied(value: float, spec: AnnotatorSpec) -> bool:
    """`_value_satisfies` against one annotator's own cutoff/direction."""
    return _value_satisfies(value, spec.cutoff, spec.direction)


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
    """Assemble one round's candidate table: smiles, source, annotator values, cutoffs.

    `cutoff_satisfied` is the primary's (NaN without a primary or a cutoff); the
    guard-rails decide `guardrails_satisfied`.
    """
    primary = _primary(annotators)
    judged = primary is not None and has_cutoff(primary)
    guards = _guards(annotators)
    rows = []
    for smi, source in source_by_smiles.items():
        row = {"smiles": smi, "source": source, **scores[smi]}
        row["cutoff_satisfied"] = (
            _cutoff_satisfied(scores[smi][primary.annotator_id], primary)
            if judged
            else float("nan")  # as an empty cell reads back
        )
        row["guardrails_satisfied"] = _cutoffs_satisfied_count(
            scores[smi], guards
        ) == len(guards)
        if tanimoto is not None:
            row["tanimoto_to_seed"] = tanimoto.get(smi, float("nan"))
        rows.append(row)
    return pd.DataFrame(rows)


def _sort_round_table(
    round_table: pd.DataFrame, primary_id: str, direction: Direction
) -> pd.DataFrame:
    """Order a round's candidates best-first by the primary.

    Ascending if its direction is "lower", else descending. NaN always sorts last
    (pandas default), so a candidate the primary couldn't score is never treated
    as "best".
    """
    if round_table.empty:
        return round_table
    return round_table.sort_values(
        primary_id, ascending=direction == "lower", kind="stable"
    ).reset_index(drop=True)


def _eligible(round_table: pd.DataFrame) -> pd.Series:
    """Candidates allowed to win: they satisfy every guard-rail (a missing score never does)."""
    return round_table["guardrails_satisfied"].astype(bool)


def _improved(
    candidate_score: float,
    best_score: float,
    direction: Direction,
    tolerance: float,
) -> bool:
    """Whether `candidate_score` beats `best_score` by more than `tolerance`."""
    if direction == "lower":
        return candidate_score < best_score - tolerance
    return candidate_score > best_score + tolerance


def _run_rounds(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    best_smiles: str | None,
    best_score: float,
    pick_winner,
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    """The loop: generate -> score -> pick_winner -> compare.

    A round whose table is already in `checkpoint` is read back, not run again.
    """
    primary = _primary(annotators)
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    columns = _table_columns(annotators, seed_smiles is not None)

    try:
        for round_num in range(1, n_rounds + 1):
            round_table = _saved_round(checkpoint, round_num, columns)
            if round_table is not None:
                logger.info(
                    f"Round {round_num}: {len(round_table)} candidates, from "
                    f"round{round_num}.csv"
                )
            else:
                round_table = _new_round(
                    generator, annotators, seed_smiles, round_num, best_smiles
                )
                _save_round(checkpoint, round_num, round_table)
            candidates_by_round[round_num] = round_table

            if round_table.empty:
                logger.warning(f"Round {round_num}: no candidates generated, stopping.")
                break

            winner = pick_winner(round_table)
            if winner is None:
                logger.warning(f"Round {round_num}: no eligible candidate, stopping.")
                break
            round_best_smiles, round_best_score = winner

            is_new_best = _improved(
                round_best_score, best_score, primary.direction, tolerance
            )
            history.append(
                {
                    "round": round_num,
                    "smiles": round_best_smiles,
                    "score": round_best_score,
                    "is_new_best": is_new_best,
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
                best_smiles, best_score = round_best_smiles, round_best_score

                # an infinite cutoff, like none (NaN), never stops a run early
                if math.isfinite(primary.cutoff) and _cutoff_satisfied(
                    best_score, primary
                ):
                    logger.info(
                        f"{primary.annotator_id}: cutoff {primary.cutoff} satisfied, "
                        "stopping."
                    )
                    break
            else:
                logger.info(f"Round {round_num}: no improvement, stopping.")
                break

    except _STOPS as error:
        _keep_finished(error, history, candidates_by_round)
        raise

    return history, candidates_by_round


def _new_round(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    round_num: int,
    best_smiles: str | None,
) -> pd.DataFrame:
    """Run one round's models: candidates around `best_smiles`, scored and sorted."""
    if logger._quiet_mode:
        logger.logger.setLevel(logging.ERROR)
    try:
        with execution.stage(f"Round {round_num} · sampling"):
            by_model = generator.generate_by_model(best_smiles)
    finally:  # a run stopped here must not leave the log silenced
        if logger._quiet_mode:
            logger.logger.setLevel(logging.INFO)
    if logger._quiet_mode:
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

    with execution.stage(f"Round {round_num} · annotation"):
        scores = _score_candidates(list(joined_source), annotators)
    tanimoto = (
        tanimoto_to_seed(seed_smiles, list(joined_source))
        if seed_smiles is not None
        else None
    )
    round_table = _round_table(joined_source, scores, annotators, tanimoto)
    primary = _primary(annotators)
    return _sort_round_table(round_table, primary.annotator_id, primary.direction)


def _run(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    best_smiles: str | None,
    best_row: dict[str, float],
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    """The rounds after the seed's: the primary optimized, the guard-rails enforced."""
    primary = _primary(annotators)
    aid = primary.annotator_id
    if best_smiles is None:
        best_score = float("-inf") if primary.direction == "higher" else float("inf")
    else:
        best_score = best_row.get(aid, float("nan"))
        if math.isnan(best_score):
            raise ValueError(f"{best_smiles!r} could not be scored by {aid!r}")

    def pick_winner(round_table: pd.DataFrame):
        eligible = round_table[_eligible(round_table) & round_table[aid].notna()]
        if eligible.empty:
            return None
        idx = (
            eligible[aid].idxmax()
            if primary.direction == "higher"
            else eligible[aid].idxmin()
        )
        row = eligible.loc[idx]
        if has_cutoff(primary) and not _cutoff_satisfied(row[aid], primary):
            if logger._quiet_mode:
                logger.warning(
                    f"{aid}: {row[aid]:.2f} (cutoff {primary.cutoff} not met)"
                )
            else:
                logger.warning(
                    f"{aid}: best eligible candidate scores {row[aid]}, which does "
                    f"not satisfy the cutoff of {primary.cutoff} (direction: "
                    f"{primary.direction})."
                )
        return row["smiles"], row[aid]

    return _run_rounds(
        generator,
        annotators,
        seed_smiles,
        n_rounds,
        tolerance,
        best_smiles=best_smiles,
        best_score=best_score,
        pick_winner=pick_winner,
        checkpoint=checkpoint,
    )
