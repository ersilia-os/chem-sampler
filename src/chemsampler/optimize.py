import dataclasses
import logging
import math
import os
from typing import Literal

import pandas as pd

from . import execution
from .models.annotator import NEEDS_SEED, TanimotoAnnotator
from .models.generator import SeedRequiredError
from .models.similarity import HubSimilarity
from .models.spec import AnnotatorDirection, AnnotatorSpec, Cutoff, Direction
from .utils.logging import logger
from .utils.similarity import tanimoto_to_seed

_RESERVED_COLUMNS = {
    "smiles",
    "source",
    "cutoffs_satisfied",
    "guardrails_satisfied",
    "weighted_score",
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


def _table_columns(
    annotators: list[AnnotatorSpec],
    tanimoto: bool,
    weighted: bool,
) -> set[str]:
    """The columns `_round_table` gives a non-empty round, with these options."""
    columns = {"smiles", "source", "cutoffs_satisfied", "guardrails_satisfied"}
    columns |= {spec.annotator_id for spec in annotators}
    if weighted:
        columns.add("weighted_score")
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
    mode: Literal["joint", "sequential", "weighted", "incremental"] | None = None,
    similarity: HubSimilarity | None = None,
    seed_smiles: str | None = None,
    n_rounds: int = 5,
    n_steps: int | None = None,
    tolerance: float = 0.0,
    checkpoint: str | os.PathLike | None = None,
) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    """
    Iteratively improve a seed molecule by alternating generation and scoring.

    Annotators are either *primary* (optimized) or *guard-rails*. A
    guard-rail is a hard constraint: a candidate that does not satisfy
    its cutoff or range is never eligible to win, in any mode, in any round or
    stage, and it takes no part in the optimization itself. The primaries are
    what `mode` combines into a round's winner:

    - "sequential": annotators are optimized one at a time, in list order. Each
      finished stage's own cutoff becomes a hard floor for every later stage -
      a later stage can trade within an already-cleared objective's accepted
      range, but never push it back below the goal declared for it. If a
      stage's round budget runs out before it ever reaches its own cutoff, its
      actual achieved value becomes the floor instead (there is nothing higher
      to promise). A stage whose own cutoff nothing can satisfy, not even the
      entering candidate, stops the whole run.
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
    - "incremental": every annotator is advanced a little at a time instead of
      all the way to its own cutoff at once. Requires `seed_smiles` and
      `n_steps`. Cycles through `annotators` in list order, `n_steps` times;
      each of these `n_steps * len(annotators)` steps is a `sequential`-style
      stage (same round budget, same early-stop-on-target, same "nothing
      clears it but something respects every floor" fallback) targeting not
      the annotator's own `cutoff` but an interim value that interpolates
      linearly between its seed score and that cutoff: `seed_value +
      (level / n_steps) * (cutoff - seed_value)`, where `level` is 1 on that
      annotator's first visit, 2 on its second, and so on - always anchored to
      the *original* seed score, never to whatever was actually achieved
      along the way. `level == n_steps` (an annotator's last visit) always
      lands exactly on its real `cutoff`. Once a step's own interim target is
      reached, that target (not the achieved value) becomes a hard floor for
      every later step, exactly as a finished stage's cutoff does in
      "sequential". Slower than "sequential" (many small stages instead of
      one per annotator), intended for difficult cases where committing an
      annotator all the way to its final cutoff before touching the next one
      risks overshooting in a way another, uncorrelated annotator can't
      recover from.

    A `TanimotoAnnotator` (`rdkit,tanimoto`) scores each candidate's similarity
    to `seed_smiles` and works like any other annotator: as a guard-rail it
    keeps candidates near the seed (or away from it), as a primary it is
    optimized.

    Parameters
    ----------
    generator : object or None
        Any object exposing `.generate_by_model(seed_smiles: str | None) ->
        dict[str, list[str]]` (e.g. `HubGenerator` or a `GeneratorPool`). May be
        `None` when `similarity` is given.
    annotators : list[AnnotatorSpec]
        At least one must be primary. In "sequential" and "incremental" modes,
        the order of the primaries is the priority/cycling order. Guard-rails:
        see above.
    mode : {"joint", "sequential", "weighted", "incremental"}, optional
        How the primaries combine: needed with several primaries; with one the
        run is "sequential" and no other mode is accepted (see `resolve_mode`).
    similarity : HubSimilarity, optional
        A similarity-search model whose nearest neighbours of the current best
        molecule join each round's candidates next to the generator's, or are
        the only candidates when `generator` is `None`. Kept separate from
        `generator`, which only takes generative models. Skipped with a warning
        when there is no seed to search from.
    seed_smiles : str, optional
        SMILES of the starting molecule. If `None`, the first round's winner is
        unconditionally the new best. Required if `mode="incremental"` or a
        `TanimotoAnnotator` is among `annotators`.
    n_rounds : int, optional
        Maximum rounds to run, by default 5. In "sequential" mode this budget
        applies separately to each stage; in "incremental" mode, separately to
        each step.
    n_steps : int, optional
        Number of interpolation levels per annotator in "incremental" mode -
        required (and must be >= 1) only for that mode; a `ValueError` if given
        with any other mode. See the "incremental" description above.
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
        One row per round attempted, with columns `round`, `smiles`, `score`,
        `is_new_best`, `active_annotator_id` (the annotator that round was
        driving; `None` for every "joint"/"weighted"-mode row), and
        `step_level` (the 1-indexed interpolation level that round belonged
        to; only present when `mode="incremental"`).
    candidates_by_round : dict[int, pandas.DataFrame]
        Every candidate considered each round (round 0 is the seed, if given),
        continuously numbered across stage/step boundaries in "sequential" and
        "incremental" modes. Columns: `smiles`, `source`, one column per
        `annotators` entry's `annotator_id`, `cutoffs_satisfied` (always
        present, regardless of `mode`; counts the primaries only),
        `guardrails_satisfied` (True when every guard-rail is
        satisfied; only such candidates can win), `weighted_score` (only present when
        `mode="weighted"`), and `tanimoto_to_seed` (only if a seed was given;
        always against `seed_smiles` itself, never the rolling best-so-far
        molecule).

    Raises
    ------
    ValueError
        If neither `generator` nor `similarity` is given, `annotators` is empty
        or has no primary one, has duplicate or reserved `annotator_id` values,
        `mode` is not one `resolve_mode` accepts, a `TanimotoAnnotator` is given
        without a seed, an entering candidate can't be scored by the annotator
        whose stage it's entering, `mode="incremental"` is given without
        `seed_smiles` or a valid `n_steps` (>= 1), `n_steps` is given with any
        other mode, or (in "incremental" mode) any annotator has a non-finite
        cutoff or can't score the seed.
    TypeError
        If `generator` is a `HubSimilarity`: pass it as `similarity` instead.

    Notes
    -----
    A run stopped midway by Ctrl-C (`KeyboardInterrupt`) or a failing model
    (`RuntimeError`, `ValueError`) re-raises that exception with the rounds
    finished before it as `partial`: the `(summary, candidates_by_round)` pair
    this function would have returned. It is absent when no round finished.

    With `checkpoint`, nothing but the round tables needs saving to continue: a
    round's winner, the floors and stages of "sequential" and "incremental",
    and when to stop are all worked out from the tables, so rounds read back
    lead to exactly the state the first call reached. The round that was
    interrupted is run again whole. Calling again with other arguments than
    the first time gives a run that neither call describes.
    """
    if isinstance(generator, HubSimilarity):
        raise TypeError(
            f"{generator.model_id} is a similarity-search model, not a generator; "
            "pass it as hill_climb(similarity=...)"
        )
    if generator is None and similarity is None:
        raise ValueError("give a generator, a similarity model, or both")
    _validate_annotators(annotators)
    mode = resolve_mode(mode, annotators)
    if mode == "incremental":
        if seed_smiles is None:
            raise ValueError("mode='incremental' requires seed_smiles")
        if n_steps is None or n_steps < 1:
            raise ValueError(
                f"mode='incremental' requires n_steps >= 1, got {n_steps!r}"
            )
    elif n_steps is not None:
        raise ValueError("n_steps is only used by mode='incremental'")
    annotators = _with_seed(annotators, seed_smiles, mode)

    if similarity is not None:
        generator = _CandidateSources(generator, similarity)

    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    best_row: dict[str, float] = {}
    best_score = float("-inf")

    try:
        if seed_smiles is not None:
            columns = _table_columns(annotators, True, mode == "weighted")
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
                    include_weighted_score=(mode == "weighted"),
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
            if mode == "incremental":
                for spec in _primaries(annotators):
                    if not math.isfinite(spec.cutoff):
                        raise ValueError(
                            f"mode='incremental' requires a finite cutoff for every "
                            f"annotator; {spec.annotator_id!r} has cutoff {spec.cutoff!r}"
                        )
                    if math.isnan(best_row[spec.annotator_id]):
                        raise ValueError(
                            f"{seed_smiles!r} could not be scored by "
                            f"{spec.annotator_id!r}, required to anchor the "
                            "mode='incremental' schedule"
                        )
            active_id = (
                None
                if mode in ("joint", "weighted")
                else _primaries(annotators)[0].annotator_id
            )
            if mode == "joint":
                best_score = seed_row["cutoffs_satisfied"]
            elif mode == "weighted":
                best_score = seed_row["weighted_score"]
            else:
                best_score = seed_row[active_id]
            round0_history = {
                "round": 0,
                "smiles": seed_smiles,
                "score": best_score,
                "is_new_best": True,
                "active_annotator_id": active_id,
            }
            if mode == "incremental":
                round0_history["step_level"] = 0
            history.append(round0_history)
            logger.info(f"Round 0: seed score = {best_score}")

        if mode == "joint":
            stage_history, stage_candidates = _run_joint(
                generator,
                annotators,
                seed_smiles,
                n_rounds,
                tolerance,
                best_smiles=seed_smiles,
                best_score=best_score,
                best_row=best_row,
                checkpoint=checkpoint,
            )
        elif mode == "weighted":
            stage_history, stage_candidates = _run_weighted(
                generator,
                annotators,
                seed_smiles,
                n_rounds,
                tolerance,
                best_smiles=seed_smiles,
                best_score=best_score,
                best_row=best_row,
                checkpoint=checkpoint,
            )
        elif mode == "incremental":
            stage_history, stage_candidates = _run_incremental(
                generator,
                annotators,
                seed_smiles,
                n_rounds,
                tolerance,
                best_smiles=seed_smiles,
                best_row=best_row,
                n_steps=n_steps,
                checkpoint=checkpoint,
            )
        else:
            stage_history, stage_candidates = _run_sequential(
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
        At least one, primary or guard-rail. See `hill_climb` for the shared
        invariants (unique, non-reserved `annotator_id` values).
    seed_smiles : str, optional
        The seed a `TanimotoAnnotator` compares `smiles` with.

    Returns
    -------
    pandas.DataFrame
        Exactly one row, shaped like a `hill_climb` round table: `smiles`,
        `source` (always "input"), one column per `annotators` entry's
        `annotator_id` (`NaN` where that annotator can't score `smiles`),
        `cutoffs_satisfied` (primaries only) and `guardrails_satisfied`.

    Raises
    ------
    ValueError
        If `annotators` is empty or has duplicate/reserved `annotator_id`
        values, or has a `TanimotoAnnotator` and no `seed_smiles` is given.
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


class _CandidateSources:
    """A generator and a similarity model behind `hill_climb`'s one-generator interface."""

    def __init__(self, generator, similarity: HubSimilarity):
        self._generator = generator
        self._similarity = similarity
        # The per-round log line reads `.generators`; a lone generator has none.
        generators = (
            [] if generator is None else getattr(generator, "generators", [generator])
        )
        self.generators = [*generators, similarity]

    def generate_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        # the generators and the similarity search run side by side on a cluster
        calls = [] if self._generator is None else [self._generator.generate_by_model]
        calls.append(self._search_by_model)
        by_model = {}
        for found in execution.current().map(lambda call: call(seed_smiles), calls):
            by_model.update(found)
        return by_model

    def _search_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        try:
            return self._similarity.search_by_model(seed_smiles)
        except SeedRequiredError:
            logger.warning(
                f"{self._similarity.model_id}: requires a seed molecule, skipped"
            )
            return {}


MODES = ("joint", "sequential", "weighted", "incremental")


def resolve_mode(mode: str | None, annotators: list[AnnotatorSpec]) -> str:
    """
    The mode a run uses: `mode` with several primaries, "sequential" with one.

    The modes differ only in how several primaries combine, so with several a
    mode must be chosen, and with one the run is "sequential", which optimizes
    that annotator until its cutoff is met.

    Parameters
    ----------
    mode : str or None
        One of `MODES`, or None. With one primary, only None or "sequential"
        (what a resumed run recorded).
    annotators : list[AnnotatorSpec]
        The run's annotators.

    Returns
    -------
    str
        One of `MODES`.

    Raises
    ------
    ValueError
        If another mode than "sequential" is given with one primary, none is
        given with several, or it is not one of `MODES`.
    """
    count = len(_primaries(annotators))
    if count == 1:
        if mode not in (None, "sequential"):
            raise ValueError(
                "with one primary annotator a run is always sequential: give no "
                "mode (--mode)"
            )
        return "sequential"
    if mode is None:
        raise ValueError(
            f"with {count} primary annotators a mode is needed (--mode): "
            + ", ".join(MODES)
        )
    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}; got {mode!r}")
    return mode


def _primaries(annotators: list[AnnotatorSpec]) -> list[AnnotatorSpec]:
    """The annotators a mode optimizes, in list order."""
    return [spec for spec in annotators if spec.priority == "primary"]


def _guards(annotators: list[AnnotatorSpec]) -> list[AnnotatorSpec]:
    """The guard-rails: constraints every candidate must satisfy."""
    return [spec for spec in annotators if spec.priority == "guard-rail"]


def _validate_annotators(
    annotators: list[AnnotatorSpec], need_primary: bool = True
) -> None:
    """Check the list-level invariants hill_climb depends on."""
    if not annotators:
        raise ValueError("annotators must not be empty")
    if need_primary and not _primaries(annotators):
        raise ValueError(
            "at least one primary annotator is required: guard-rails are not optimized"
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
    annotators: list[AnnotatorSpec], seed_smiles: str | None, mode: str | None = None
) -> list[AnnotatorSpec]:
    """The annotators, with every Tanimoto one comparing candidates to the seed."""
    tanimoto = [s for s in annotators if isinstance(s.annotator, TanimotoAnnotator)]
    if not tanimoto:
        return annotators
    if seed_smiles is None:
        raise ValueError(NEEDS_SEED)
    primaries = _primaries(annotators)
    first = primaries[0] if primaries else None
    if (
        mode in ("sequential", "incremental")
        and first in tanimoto
        and first.direction == "higher"
    ):
        logger.warning(
            f"{first.annotator_id} (higher) is the first primary: nothing can be "
            "more similar to the seed than the seed itself (1.0), so its stage "
            "cannot improve; put it after another primary"
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


def _cutoff_satisfied_mask(
    values: pd.Series, cutoff: Cutoff, direction: AnnotatorDirection
) -> pd.Series:
    """Vectorized `_cutoff_satisfied`; NaN comparisons are already False in pandas."""
    if direction == "range":
        low, high = cutoff
        return (values >= low) & (values <= high)
    if direction == "higher":
        return values >= cutoff
    return values <= cutoff


def _cutoffs_satisfied_count(
    scores_row: dict[str, float], annotators: list[AnnotatorSpec]
) -> int:
    return sum(
        _cutoff_satisfied(scores_row[spec.annotator_id], spec) for spec in annotators
    )


def _incremental_target(
    seed_value: float, cutoff: float, level: int, n_steps: int
) -> float:
    """Linear interpolation target for one annotator at one step level (1..n_steps).

    Direction-agnostic: works whether `cutoff` is above or below `seed_value`.
    `level == n_steps` always lands exactly on `cutoff`.
    """
    return seed_value + (level / n_steps) * (cutoff - seed_value)


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
    include_weighted_score: bool = False,
) -> pd.DataFrame:
    """Assemble one round's candidate table: smiles, source, annotator values, cutoff count.

    The count and the weighted score cover the primary annotators only; the
    guard-rails decide `guardrails_satisfied`.
    """
    primaries = _primaries(annotators)
    guards = _guards(annotators)
    rows = []
    for smi, source in source_by_smiles.items():
        row = {"smiles": smi, "source": source, **scores[smi]}
        row["cutoffs_satisfied"] = _cutoffs_satisfied_count(scores[smi], primaries)
        row["guardrails_satisfied"] = _cutoffs_satisfied_count(
            scores[smi], guards
        ) == len(guards)
        if include_weighted_score:
            row["weighted_score"] = _weighted_score(scores[smi], primaries)
        if tanimoto is not None:
            row["tanimoto_to_seed"] = tanimoto.get(smi, float("nan"))
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


def _eligible(round_table: pd.DataFrame) -> pd.Series:
    """Candidates allowed to win: they satisfy every guard-rail (a missing score never does)."""
    return round_table["guardrails_satisfied"].astype(bool)


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
    active_cutoff: float | None = None,
    sort_column: str = "cutoffs_satisfied",
    include_weighted_score: bool = False,
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame], str | None, float, dict[str, float]]:
    """Shared loop: generate -> score -> pick_winner -> compare. Used by every mode.

    A round whose table is already in `checkpoint` is read back, not run again.
    """
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    columns = _table_columns(
        annotators,
        seed_smiles is not None,
        include_weighted_score,
    )

    try:
        for offset in range(n_rounds):
            round_num = start_round + offset

            round_table = _saved_round(checkpoint, round_num, columns)
            if round_table is not None:
                logger.info(
                    f"Round {round_num}: {len(round_table)} candidates, from "
                    f"round{round_num}.csv"
                )
            else:
                round_table = _new_round(
                    generator,
                    annotators,
                    seed_smiles,
                    round_num,
                    best_smiles,
                    active_annotator_id,
                    winner_direction,
                    sort_column,
                    include_weighted_score,
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

    except _STOPS as error:
        _keep_finished(error, history, candidates_by_round)
        raise

    return history, candidates_by_round, best_smiles, best_score, best_row


def _new_round(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    round_num: int,
    best_smiles: str | None,
    active_annotator_id: str | None,
    winner_direction: Direction | None,
    sort_column: str,
    include_weighted_score: bool,
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
    round_table = _round_table(
        joined_source,
        scores,
        annotators,
        tanimoto,
        include_weighted_score,
    )
    round_table = _sort_round_table(
        round_table, active_annotator_id, winner_direction, sort_column
    )
    return round_table


def _row_values(row: pd.Series, annotators: list[AnnotatorSpec]) -> dict[str, float]:
    """Every annotator's value for one round-table row, as a plain dict."""
    return {spec.annotator_id: row[spec.annotator_id] for spec in annotators}


def _run_joint(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    def pick_winner(round_table: pd.DataFrame):
        eligible = round_table[_eligible(round_table)]
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
        checkpoint=checkpoint,
    )
    return history, candidates_by_round


def _run_weighted(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    best_smiles: str | None,
    best_score: float,
    best_row: dict[str, float],
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    def pick_winner(round_table: pd.DataFrame):
        mask = _eligible(round_table)
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
        include_weighted_score=True,
        checkpoint=checkpoint,
    )
    return history, candidates_by_round


def _run_sequential(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    best_smiles: str | None,
    best_row: dict[str, float],
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    floors: list[tuple[str, float, Direction]] = []
    next_round = 1

    try:
        for spec in _primaries(annotators):
            if best_smiles is None:
                best_score = (
                    float("-inf") if spec.direction == "higher" else float("inf")
                )
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
                mask = _eligible(round_table)
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
                    floor_mask = _eligible(round_table)
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
                return (
                    row["smiles"],
                    row[spec.annotator_id],
                    _row_values(row, annotators),
                )

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
                active_cutoff=spec.cutoff,
                checkpoint=checkpoint,
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

            floor_value = spec.cutoff if _cutoff_satisfied(achieved, spec) else achieved
            floors.append((spec.annotator_id, floor_value, spec.direction))

    except _STOPS as error:
        _keep_finished(error, history, candidates_by_round)
        raise

    return history, candidates_by_round


def _run_incremental(
    generator,
    annotators: list[AnnotatorSpec],
    seed_smiles: str,
    n_rounds: int,
    tolerance: float,
    best_smiles: str,
    best_row: dict[str, float],
    n_steps: int,
    checkpoint=None,
) -> tuple[list[dict], dict[int, pd.DataFrame]]:
    """Like `_run_sequential`, but cycling `annotators` `n_steps` times through
    interpolated interim targets instead of once through each own cutoff."""
    history = []
    candidates_by_round: dict[int, pd.DataFrame] = {}
    seed_value = dict(best_row)
    floors: dict[str, tuple[float, Direction]] = {}
    next_round = 1

    level = None
    try:
        for level in range(1, n_steps + 1):
            for spec in _primaries(annotators):
                target = _incremental_target(
                    seed_value[spec.annotator_id], spec.cutoff, level, n_steps
                )

                entry_value = best_row.get(spec.annotator_id, float("nan"))
                if math.isnan(entry_value):
                    raise ValueError(
                        f"{best_smiles!r} could not be scored by {spec.annotator_id!r} "
                        "entering its step"
                    )
                best_score = entry_value
                entry_ok = _value_satisfies(entry_value, target, spec.direction)

                floors_snapshot = [
                    (aid, floor_value, floor_direction)
                    for aid, (floor_value, floor_direction) in floors.items()
                    if aid != spec.annotator_id
                ]

                def pick_winner(
                    round_table: pd.DataFrame,
                    spec=spec,
                    target=target,
                    floors_snapshot=floors_snapshot,
                ):
                    mask = _eligible(round_table)
                    mask &= _cutoff_satisfied_mask(
                        round_table[spec.annotator_id], target, spec.direction
                    )
                    for floor_id, floor_value, floor_direction in floors_snapshot:
                        mask &= _cutoff_satisfied_mask(
                            round_table[floor_id], floor_value, floor_direction
                        )
                    eligible = round_table[mask]
                    if eligible.empty:
                        # No candidate clears this step's own target even with floors
                        # intact. Still report an improvement if one exists that
                        # respects every prior floor - floors are a hard invariant and
                        # are never dropped, even in this fallback (own target is
                        # still excluded; that's the whole point of the fallback).
                        floor_mask = _eligible(round_table)
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
                                    f"{spec.annotator_id}: {row[spec.annotator_id]:.2f} "
                                    f"(step target {target} not met)"
                                )
                            else:
                                logger.warning(
                                    f"{spec.annotator_id}: improved to "
                                    f"{row[spec.annotator_id]}, but did not satisfy "
                                    f"this step's target of {target} (direction: "
                                    f"{spec.direction}). Continuing to next round."
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
                    return (
                        row["smiles"],
                        row[spec.annotator_id],
                        _row_values(row, annotators),
                    )

                (
                    step_history,
                    step_candidates,
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
                    active_cutoff=target,
                    checkpoint=checkpoint,
                )
                for row in step_history:
                    row["step_level"] = level
                history.extend(step_history)
                candidates_by_round.update(step_candidates)
                next_round += len(step_candidates)

                if new_best_smiles is None:
                    logger.warning(
                        f"{spec.annotator_id}: no candidate found, stopping the whole run."
                    )
                    return history, candidates_by_round

                if new_best_smiles == best_smiles:
                    if not entry_ok:
                        logger.warning(
                            f"{spec.annotator_id}: no candidate satisfies this step's "
                            "target, stopping the whole run."
                        )
                        return history, candidates_by_round
                    achieved = best_score
                else:
                    achieved = new_best_score
                    best_smiles, best_row = new_best_smiles, new_best_row

                floor_value = (
                    target
                    if _value_satisfies(achieved, target, spec.direction)
                    else achieved
                )
                floors[spec.annotator_id] = (floor_value, spec.direction)

    except _STOPS as error:
        if level is not None:  # rounds of the interrupted step lack their level
            for row in getattr(error, "_finished_rounds", ([], {}))[0]:
                row.setdefault("step_level", level)
        _keep_finished(error, history, candidates_by_round)
        raise

    return history, candidates_by_round
