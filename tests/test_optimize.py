import math

import pandas as pd
import pytest

from chemsampler.models.annotator import TanimotoAnnotator
from chemsampler.models.generator import SeedRequiredError
from chemsampler.models.similarity import HubSimilarity
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import annotate, hill_climb, resolve_mode, write_results
from chemsampler.utils.logging import logger


class StubGenerator:
    """Returns one fixed set of candidates per round, regardless of the seed."""

    model_id = "stub"

    def __init__(self, candidates_by_round):
        self._candidates_by_round = list(candidates_by_round)
        self._call_count = 0

    def generate_by_model(self, seed_smiles):
        candidates = self._candidates_by_round[self._call_count]
        self._call_count += 1
        return {self.model_id: list(candidates)}


class StubAnnotator:
    """Scores SMILES by a fixed lookup table, dropping SMILES it has no entry for."""

    def __init__(self, scores):
        self._scores = scores

    def score(self, smiles_list):
        return {smi: self._scores[smi] for smi in smiles_list if smi in self._scores}


# --- translated from the pre-redesign suite, as mode="sequential" with a --------
# --- permissive cutoff: this is the only faithful single-annotator migration ---


def test_stops_when_a_round_does_not_improve():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 0.5}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=5
    )

    assert list(summary["round"]) == [0, 1, 2]
    assert list(summary["is_new_best"]) == [True, True, False]
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert set(candidates_by_round) == {0, 1, 2}


def test_respects_n_rounds_cap():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=2
    )

    assert list(summary["round"]) == [0, 1, 2]
    assert summary["is_new_best"].all()


def test_summary_has_expected_columns():
    generator = StubGenerator([["CCO"]])
    annotators = [
        AnnotatorSpec(
            "score",
            StubAnnotator({"CCO": 1.0}),
            cutoff=float("-inf"),
            direction="higher",
        )
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert list(summary.columns) == [
        "round",
        "smiles",
        "score",
        "is_new_best",
        "active_annotator_id",
    ]


def test_candidates_table_has_expected_columns():
    generator = StubGenerator([["CCO"]])
    annotators = [
        AnnotatorSpec(
            "score",
            StubAnnotator({"CCO": 1.0}),
            cutoff=float("-inf"),
            direction="higher",
        )
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1].columns) == [
        "smiles",
        "source",
        "score",
        "cutoffs_satisfied",
        "guardrails_satisfied",
        "tanimoto_to_seed",
    ]


def test_best_score_trajectory_is_non_decreasing():
    scores = {"CCO": 1.0, "CCC": 1.5, "CCN": 1.2, "CC": 2.0}
    generator = StubGenerator([["CCC"], ["CCN"], ["CC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=5
    )

    best_so_far = summary[summary["is_new_best"]]["score"]
    assert list(best_so_far) == sorted(best_so_far)


def test_tolerance_requires_margin_to_count_as_improvement():
    scores = {"CCO": 1.0, "CCC": 1.05}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="sequential",
        seed_smiles="CCO",
        n_rounds=5,
        tolerance=0.1,
    )

    assert list(summary["round"]) == [0, 1]
    assert list(summary["is_new_best"]) == [True, False]


def test_seedless_first_round_is_unconditionally_the_new_best():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "score",
            StubAnnotator({"CCC": -5.0}),
            cutoff=float("-inf"),
            direction="higher",
        )
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", n_rounds=1
    )

    assert list(summary["round"]) == [1]
    assert summary.iloc[0]["is_new_best"] == True
    assert 0 not in candidates_by_round
    assert "tanimoto_to_seed" not in candidates_by_round[1].columns


def test_a_tanimoto_annotator_requires_a_seed():
    annotators = [
        AnnotatorSpec("score", StubAnnotator({"CCC": 1.0}), 0.0, "higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), 0.2, "higher", priority="guard-rail"
        ),
    ]

    # refused before anything is generated: StubGenerator([]) fails if called
    with pytest.raises(ValueError, match="tanimoto compares candidates with the seed"):
        hill_climb(StubGenerator([]), annotators, seed_smiles=None)


def test_sequential_mode_single_annotator_reproduces_pure_maximization():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=2
    )

    assert list(summary["is_new_best"]) == [True, True, True]
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert summary.iloc[-1]["score"] == 3.0


# --- validation ------------------------------------------------------------


def test_one_primary_defaults_to_sequential_and_several_need_a_mode():
    scores = {"CCO": 0.0, "CCC": 1.0}
    one = [AnnotatorSpec("a", StubAnnotator(scores), float("inf"), "higher")]
    two = [*one, AnnotatorSpec("b", StubAnnotator(scores), 0.0, "higher")]

    summary, _ = hill_climb(
        StubGenerator([["CCC"]]), one, seed_smiles="CCO", n_rounds=1
    )

    assert summary["active_annotator_id"].dropna().unique().tolist() == ["a"]
    assert resolve_mode(None, one) == "sequential"
    assert resolve_mode("joint", one) == "joint"
    with pytest.raises(ValueError, match="with 2 primary annotators a mode is needed"):
        hill_climb(StubGenerator([]), two, seed_smiles="CCO")
    with pytest.raises(ValueError, match="mode must be one of"):
        resolve_mode("greedy", two)


def test_hill_climb_rejects_empty_annotators():
    with pytest.raises(ValueError, match="empty"):
        hill_climb(StubGenerator([]), [], mode="joint", seed_smiles="CCO")


def test_hill_climb_rejects_reserved_annotator_id():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "cutoffs_satisfied",
            StubAnnotator({"CCO": 1.0}),
            cutoff=0.0,
            direction="higher",
        )
    ]

    with pytest.raises(ValueError, match="reserved"):
        hill_climb(generator, annotators, mode="joint", seed_smiles="CCO")


# --- sequential mode: floors, own-cutoff gating, abort vs. continue --------


def test_sequential_mode_floor_is_the_cutoff_not_the_achieved_value():
    a_scores = {"CCO": 1.0, "CCC": 10.0, "CCN": 2.0, "CC": 12.0}
    b_scores = {"CCO": 0.0, "CCC": 3.0, "CCN": 8.0, "CC": 6.0}
    generator = StubGenerator([["CCC"], ["CCN", "CC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # "a" achieves 10.0 in stage 1, well above its cutoff of 0.0 - but the
    # floor for stage 2 is the *cutoff* (0.0), not the incidental 10.0 it
    # happened to land on. Both CCN (2.0) and CC (12.0) clear that floor, so
    # stage 2 is free to pick on "b" alone: CCN wins with the higher score
    # (8.0 > 6.0), even though it would have failed a floor of 10.0.
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_sequential_mode_does_not_enforce_not_yet_reached_annotators():
    a_scores = {"CCO": 1.0, "CCC": 10.0}
    b_scores = {
        "CCO": 0.0,
        "CCC": -100.0,
    }  # CCC badly fails "b", which hasn't started yet
    generator = StubGenerator([["CCC"], []])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert summary[summary["active_annotator_id"] == "a"].iloc[-1]["smiles"] == "CCC"


def test_sequential_mode_own_cutoff_gates_its_own_stage():
    a_scores = {"CCO": 1.0, "CCC": 2.0}  # CCC beats the seed but not the cutoff
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher")
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # With new behavior: improvement (1.0→2.0) is recorded even if cutoff (5.0) not met
    assert list(summary["round"]) == [0, 1]
    assert summary.iloc[1]["score"] == 2.0  # Round 1 improved
    assert (
        candidates_by_round[1].iloc[0]["cutoffs_satisfied"] == 0
    )  # But didn't satisfy cutoff


def test_sequential_stage_with_no_eligible_candidate_but_improvement_continues():
    # Test that improvement without meeting cutoff still advances to next stage
    a_scores = {"CCO": 1.0, "CCC": 2.0}
    b_scores = {"CCO": 1.0, "CCC": 3.0}
    # First round generates CCC for stage a, second round repeats (no improvement, done)
    generator = StubGenerator([["CCC"], ["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=0.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # Stage "a" improves (1.0 → 2.0) but doesn't clear its cutoff (5.0).
    # With its n_rounds budget spent, stage "b" starts and improves further (2.0 → 3.0).
    assert list(summary["active_annotator_id"]) == ["a", "a", "b"]
    assert summary.iloc[0]["score"] == 1.0  # Seed
    assert summary.iloc[1]["score"] == 2.0  # Stage a improved from seed
    assert summary.iloc[2]["score"] == 3.0  # Stage b improved further


def test_sequential_fallback_respects_prior_stage_floor():
    # "a" achieves 10.0 in stage 1, clearing its cutoff of 8.0 - so the floor
    # for stage 2 is 8.0. In stage 2, no candidate clears "b"'s (unreachable)
    # cutoff, so the no-eligible-candidate fallback kicks in. CCN has by far
    # the best "b" value but violates the floor (2.0 < 8.0); CC has a worse
    # "b" value but respects it (12.0 >= 8.0). The fallback must never trade
    # away the floor, so CC must win, not CCN.
    a_scores = {"CCO": 1.0, "CCC": 10.0, "CCN": 2.0, "CC": 12.0}
    b_scores = {"CCO": 1.0, "CCC": 5.0, "CCN": 50.0, "CC": 20.0}
    generator = StubGenerator([["CCC"], ["CCN", "CC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=8.0, direction="higher"),
        AnnotatorSpec(
            "b", StubAnnotator(b_scores), cutoff=100.0, direction="higher"
        ),  # unreachable, so stage 2 always falls back
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert summary.iloc[-1]["smiles"] == "CC"
    assert summary.iloc[-1]["score"] == 20.0


def test_sequential_fallback_stops_run_when_nothing_respects_the_floor():
    # Same shape as above (floor = 8.0, the cutoff "a" clears at 10.0), but
    # both stage-2 candidates violate it - there is no floor-respecting
    # candidate to fall back to at all, so the run must stop rather than
    # return a floor-violating candidate.
    a_scores = {"CCO": 1.0, "CCC": 10.0, "CCN": 2.0, "CC": 3.0}
    b_scores = {"CCO": 1.0, "CCC": 5.0, "CCN": 50.0, "CC": 20.0}
    generator = StubGenerator([["CCC"], ["CCN", "CC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=8.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=100.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert list(summary["active_annotator_id"]) == ["a", "a"]
    assert summary.iloc[-1]["smiles"] == "CCC"


def test_sequential_stage_with_no_improvement_but_entry_eligible_moves_to_next_stage():
    a_scores = {"CCO": 10.0, "CCC": 3.0, "CCN": 15.0}  # CCC is worse than the seed
    b_scores = {"CCO": 0.0, "CCC": 0.0, "CCN": 7.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # Stage "a" doesn't improve on the seed, but the seed already clears "a"'s
    # cutoff (10.0 >= 0.0), so the cutoff (0.0) locks in as the floor and
    # stage "b" still runs. (round 0 = seed, round 1 = "a"'s failed attempt,
    # round 2 = "b".)
    assert list(summary["active_annotator_id"]) == ["a", "a", "b"]
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_sequential_mode_candidate_missing_active_score_cannot_win():
    # "CCC" is deliberately absent from "a"'s table.
    a_scores = {"CCO": 1.0, "CCN": 0.5}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec(
            "a", StubAnnotator(a_scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert math.isnan(round1.loc["CCC", "a"])
    assert summary.iloc[-1]["smiles"] == "CCN"


# --- incremental mode: interpolated targets, floor-per-level, validation ---


def test_incremental_mode_requires_seed_smiles():
    annotators = [
        AnnotatorSpec("a", StubAnnotator({"CCO": 1.0}), cutoff=10.0, direction="higher")
    ]
    with pytest.raises(ValueError, match="seed_smiles"):
        hill_climb(StubGenerator([]), annotators, mode="incremental", n_steps=5)


def test_incremental_mode_requires_n_steps():
    annotators = [
        AnnotatorSpec("a", StubAnnotator({"CCO": 1.0}), cutoff=10.0, direction="higher")
    ]
    with pytest.raises(ValueError, match="n_steps"):
        hill_climb(StubGenerator([]), annotators, mode="incremental", seed_smiles="CCO")

    with pytest.raises(ValueError, match="n_steps"):
        hill_climb(
            StubGenerator([]),
            annotators,
            mode="incremental",
            seed_smiles="CCO",
            n_steps=0,
        )


def test_n_steps_rejected_outside_incremental_mode():
    annotators = [
        AnnotatorSpec("a", StubAnnotator({"CCO": 1.0}), cutoff=10.0, direction="higher")
    ]
    with pytest.raises(ValueError, match="n_steps"):
        hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
            mode="sequential",
            seed_smiles="CCO",
            n_steps=3,
        )


def test_incremental_mode_requires_finite_cutoff():
    annotators = [
        AnnotatorSpec(
            "a", StubAnnotator({"CCO": 1.0}), cutoff=float("inf"), direction="higher"
        )
    ]
    with pytest.raises(ValueError, match="finite"):
        hill_climb(
            StubGenerator([]),
            annotators,
            mode="incremental",
            seed_smiles="CCO",
            n_steps=5,
        )


def test_incremental_mode_requires_seed_scoreable_by_every_annotator():
    annotators = [
        AnnotatorSpec("a", StubAnnotator({}), cutoff=10.0, direction="higher")
    ]  # "CCO" missing from the score table
    with pytest.raises(ValueError, match="anchor"):
        hill_climb(
            StubGenerator([]),
            annotators,
            mode="incremental",
            seed_smiles="CCO",
            n_steps=5,
        )


def test_incremental_mode_interpolates_and_cycles_annotators_round_robin():
    # seed maip=10 -> cutoff 50 (higher); seed mw=100 -> cutoff 60 (lower);
    # n_steps=2. Level-1 targets: maip 10+(1/2)*(50-10)=30, mw
    # 100+(1/2)*(60-100)=80. Level-2 targets equal the real cutoffs (50, 60).
    maip_scores = {
        "CCO": 10.0,
        "A1": 29.9,
        "A2": 30.0,  # round 1 (maip, level 1, target 30): only A2 clears it
        "A3": 40.0,
        "A4": 35.0,  # round 2: floor check only (need >= 30)
        "A5": 49.9,
        "A6": 50.0,  # round 3 (maip, level 2, target 50): only A6 clears it
        "A7": 55.0,
        "A8": 55.0,  # round 4: floor check only (need >= 50)
    }
    mw_scores = {
        "CCO": 100.0,
        "A1": 90.0,
        "A2": 85.0,  # round 1: not gated yet, values irrelevant
        "A3": 80.1,
        "A4": 80.0,  # round 2 (mw, level 1, target 80): only A4 clears it
        "A5": 75.0,
        "A6": 75.0,  # round 3: floor check only (need <= 80)
        "A7": 60.1,
        "A8": 60.0,  # round 4 (mw, level 2, target 60): only A8 clears it
    }
    generator = StubGenerator([["A1", "A2"], ["A3", "A4"], ["A5", "A6"], ["A7", "A8"]])
    annotators = [
        AnnotatorSpec(
            "maip", StubAnnotator(maip_scores), cutoff=50.0, direction="higher"
        ),
        AnnotatorSpec("mw", StubAnnotator(mw_scores), cutoff=60.0, direction="lower"),
    ]

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="incremental",
        seed_smiles="CCO",
        n_rounds=1,
        n_steps=2,
    )

    assert list(summary["smiles"]) == ["CCO", "A2", "A4", "A6", "A8"]
    assert list(summary["active_annotator_id"]) == ["maip", "maip", "mw", "maip", "mw"]
    assert list(summary["step_level"]) == [0, 1, 1, 2, 2]
    assert summary.iloc[-1]["score"] == 60.0  # final level lands exactly on mw's cutoff


def test_incremental_mode_floor_is_the_target_not_the_achieved_value():
    # n_steps=1, so each annotator's single visit targets its real cutoff
    # (maip >= 50.0, mw <= 60.0). "maip" overshoots to 100.0 - but the floor
    # for "mw"'s step must be the *target* (50.0), not the incidental 100.0.
    # "Y1" (maip=70.0, mw=55.0) clears a floor of 50.0 and improves on X1's mw
    # (95.0), but would fail a floor of 100.0, so it only wins if the floor
    # was computed correctly.
    maip_scores = {"CCO": 10.0, "X1": 100.0, "X2": 20.0, "Y1": 70.0}
    mw_scores = {"CCO": 100.0, "X1": 95.0, "Y1": 55.0}
    generator = StubGenerator([["X1", "X2"], ["Y1"]])
    annotators = [
        AnnotatorSpec(
            "maip", StubAnnotator(maip_scores), cutoff=50.0, direction="higher"
        ),
        AnnotatorSpec("mw", StubAnnotator(mw_scores), cutoff=60.0, direction="lower"),
    ]

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="incremental",
        seed_smiles="CCO",
        n_rounds=1,
        n_steps=1,
    )

    assert list(summary["smiles"]) == ["CCO", "X1", "Y1"]


def test_incremental_mode_stops_whole_run_when_a_step_finds_no_candidates():
    annotators = [
        AnnotatorSpec("a", StubAnnotator({"CCO": 1.0}), cutoff=10.0, direction="higher")
    ]
    generator = StubGenerator([[]])

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="incremental",
        seed_smiles="CCO",
        n_rounds=1,
        n_steps=3,
    )

    assert list(summary["smiles"]) == ["CCO"]


def test_incremental_mode_stops_whole_run_when_step_has_no_improvement_and_entry_unmet():
    a_scores = {"CCO": 1.0, "CCC": 0.5}  # CCC is worse than the seed
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=10.0, direction="higher")
    ]

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="incremental",
        seed_smiles="CCO",
        n_rounds=1,
        n_steps=3,
    )

    # The failed attempt ("CCC", worse than the seed) is still recorded as an
    # attempted round, but is_new_best=False and the run stops there - it
    # never reaches level 2 or 3 (which would need a 2nd/3rd generator round).
    assert list(summary["smiles"]) == ["CCO", "CCC"]
    assert list(summary["is_new_best"]) == [True, False]
    assert list(summary["step_level"]) == [0, 1]


def test_step_level_column_only_present_in_incremental_mode():
    scores = {"CCO": 1.0, "CCC": 2.0}
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher")
    ]

    for mode in ("sequential", "joint", "weighted"):
        summary, _ = hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
            mode=mode,
            seed_smiles="CCO",
            n_rounds=1,
        )
        assert "step_level" not in summary.columns

    summary, _ = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        mode="incremental",
        seed_smiles="CCO",
        n_rounds=1,
        n_steps=1,
    )
    assert "step_level" in summary.columns


# --- joint mode: count-based scoring, tie-breaking, tanimoto exclusion ----


def test_joint_mode_single_annotator_is_threshold_not_maximization():
    scores = {"CCO": 1.0, "CCC": 10.0, "CCN": 6.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(scores), cutoff=5.0, direction="higher")
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    # Both CCC and CCN clear the single cutoff, so joint mode can't tell them
    # apart by raw score (both score cutoffs_satisfied=1) - unlike sequential
    # mode, the higher raw value (CCC, 10.0) is not guaranteed to win.
    round1 = candidates_by_round[1]
    assert set(round1["cutoffs_satisfied"]) == {1}


def test_joint_mode_round_score_is_cutoffs_satisfied_count():
    a_scores = {"CCO": 0.0, "CCC": 10.0, "CCN": 10.0}
    b_scores = {"CCO": 0.0, "CCC": 0.0, "CCN": 10.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["CCC", "cutoffs_satisfied"] == 1
    assert round1.loc["CCN", "cutoffs_satisfied"] == 2
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert summary.iloc[-1]["score"] == 2


def test_joint_mode_breaks_ties_by_first_found_candidate():
    scores = {"CCO": 0.0, "CCC": 10.0, "CCN": 10.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=5.0, direction="higher")
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    # CCC and CCN are tied (both cutoffs_satisfied=1); CCC wins as the first
    # one encountered, per the generator's own order.
    assert summary.iloc[-1]["smiles"] == "CCC"


def test_joint_mode_candidate_missing_one_annotator_only_zeroes_that_contribution():
    a_scores = {"CCO": 1.0, "CCN": 0.5}  # "CCC" missing from "a"
    b_scores = {"CCO": 1.0, "CCC": 10.0, "CCN": 10.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert math.isnan(round1.loc["CCC", "a"])
    assert round1.loc["CCC", "cutoffs_satisfied"] == 1
    assert round1.loc["CCN", "cutoffs_satisfied"] == 2
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_joint_mode_stops_when_nothing_clears_any_cutoff():
    a_scores = {"CCO": 1.0, "CCC": 0.5}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher")
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=5
    )

    assert list(summary["round"]) == [0, 1]
    assert list(summary["is_new_best"]) == [True, False]
    assert len(candidates_by_round[1]) == 1
    assert candidates_by_round[1].iloc[0]["cutoffs_satisfied"] == 0


def test_active_annotator_id_is_none_for_joint_mode():
    scores = {"CCO": 1.0, "CCC": 10.0}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher")
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    assert summary["active_annotator_id"].isna().all()


def test_a_tanimoto_guard_rail_is_not_counted_in_joint_mode():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    a_scores = {"CCO": 1.0, "CO": 1.0, "c1ccc2ccccc2c1": 1.0}
    b_scores = {"CCO": 1.0, "CO": 0.0, "c1ccc2ccccc2c1": 10.0}
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), 0.2, "higher", priority="guard-rail"
        ),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    # Naphthalene clears both cutoffs (count=2) but fails the Tanimoto guard-rail
    # (0.0 < 0.2); methanol clears only "a" (count=1) but passes it
    # (0.2857 >= 0.2) and wins despite the lower count.
    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["c1ccc2ccccc2c1", "cutoffs_satisfied"] == 2
    assert round1.loc["CO", "tanimoto"] == 0.2857142857142857
    assert candidates_by_round[0].iloc[0]["tanimoto"] == 1.0
    assert summary.iloc[-1]["smiles"] == "CO"
    assert summary.iloc[-1]["score"] == 1


# --- weighted mode: sign correction, weighting, missing values -------------


def test_weighted_mode_single_annotator_reproduces_pure_maximization():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher")
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    # Weight is irrelevant with a single annotator - the highest raw value wins.
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert summary.iloc[-1]["score"] == 3.0


def test_weighted_mode_respects_lower_direction():
    scores = {"CCO": 10.0, "CCC": 8.0, "CCN": 2.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="lower")
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    # "lower" negates the value before summing, so the *smallest* raw value
    # (CCN, 2.0) has the largest weighted_score and wins.
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_weighted_mode_combines_multiple_annotators_by_weight():
    a_scores = {"CCO": 0.0, "CCC": 10.0, "CCN": 1.0}
    b_scores = {"CCO": 0.0, "CCC": 1.0, "CCN": 10.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec(
            "a", StubAnnotator(a_scores), cutoff=0.0, direction="higher", weight=5.0
        ),
        AnnotatorSpec(
            "b", StubAnnotator(b_scores), cutoff=0.0, direction="higher", weight=1.0
        ),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    # CCC: 5*10 + 1*1 = 51. CCN: 5*1 + 1*10 = 15. CCC wins because "a" is
    # weighted 5x heavier, even though CCN has the higher "b" value and would
    # win a joint-mode (equal-weight, cutoffs-satisfied) tie-break.
    assert summary.iloc[-1]["smiles"] == "CCC"
    assert summary.iloc[-1]["score"] == 51.0


def test_weighted_mode_candidate_missing_one_annotator_is_ineligible():
    a_scores = {"CCO": 1.0, "CCN": 0.5}  # "CCC" missing from "a"
    b_scores = {"CCO": 1.0, "CCC": 100.0, "CCN": 1.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=0.0, direction="higher"),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    # CCC has by far the higher raw "b" value, but its missing "a" score makes
    # its weighted_score NaN - ineligible - so CCN wins despite the lower total.
    round1 = candidates_by_round[1].set_index("smiles")
    assert math.isnan(round1.loc["CCC", "weighted_score"])
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_active_annotator_id_is_none_for_weighted_mode():
    scores = {"CCO": 1.0, "CCC": 10.0}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher")
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    assert summary["active_annotator_id"].isna().all()


def test_a_tanimoto_guard_rail_is_not_weighted():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    scores = {"CCO": 1.0, "CO": 1.0, "c1ccc2ccccc2c1": 100.0}
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), 0.2, "higher", priority="guard-rail"
        ),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    # Naphthalene has the far higher weighted_score but fails the Tanimoto
    # guard-rail (0.0 < 0.2); methanol passes (0.2857 >= 0.2) and wins.
    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["c1ccc2ccccc2c1", "weighted_score"] == 100.0
    assert summary.iloc[-1]["smiles"] == "CO"


def test_weighted_mode_round_table_sorted_by_weighted_score_descending():
    a_scores = {"CCO": 0.0, "CCC": 1.0, "CCN": 5.0, "CC": 2.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher")
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="weighted", seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["weighted_score"]) == [5.0, 2.0, 1.0]


def test_weighted_score_column_only_present_in_weighted_mode():
    scores = {"CCO": 1.0, "CCC": 2.0}
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), cutoff=0.0, direction="higher")
    ]

    for mode in ("sequential", "joint"):
        _, candidates_by_round = hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
            mode=mode,
            seed_smiles="CCO",
            n_rounds=1,
        )
        assert "weighted_score" not in candidates_by_round[1].columns

    _, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        mode="weighted",
        seed_smiles="CCO",
        n_rounds=1,
    )
    assert "weighted_score" in candidates_by_round[1].columns


def test_a_lower_tanimoto_guard_rail_seeks_novelty():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    scores = {"CCO": 1.0, "CO": 1.0, "c1ccc2ccccc2c1": 1.0}
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), 0.2, "lower", priority="guard-rail"
        ),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    # Both candidates clear the annotator cutoff equally, so the guard-rail
    # decides: methanol (0.2857) is too similar, naphthalene (0.0) wins.
    assert summary.iloc[-1]["smiles"] == "c1ccc2ccccc2c1"


def test_a_tanimoto_range_guard_rail_keeps_similar_but_not_identical():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    scores = {"CCO": 1.0, "CCO.O": 5.0, "CO": 2.0, "c1ccc2ccccc2c1": 9.0}
    generator = StubGenerator([["CCO.O", "CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(scores), float("-inf"), "higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), (0.2, 0.8), "range", priority="guard-rail"
        ),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # "CCO.O" (ethanol and water) is too close to the seed (0.857), naphthalene
    # too far (0.0); methanol (0.286) is in between and wins.
    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["CCO.O", "tanimoto"] == 0.8571428571428571
    assert summary.iloc[-1]["smiles"] == "CO"


def test_a_lower_tanimoto_primary_drives_the_search_away_from_the_seed():
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("tanimoto", TanimotoAnnotator(), 0.1, "lower"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert summary.iloc[-1]["smiles"] == "c1ccc2ccccc2c1"
    assert summary.iloc[-1]["score"] == 0.0


def test_a_higher_tanimoto_first_primary_warns(caplog, verbose_logger):
    annotators = [
        AnnotatorSpec("tanimoto", TanimotoAnnotator(), 0.5, "higher"),
        AnnotatorSpec("score", StubAnnotator({"CCO": 1.0, "CO": 2.0}), 0.0, "higher"),
    ]

    with caplog.at_level("WARNING"):
        hill_climb(
            StubGenerator([["CO"], ["CO"]]),
            annotators,
            mode="sequential",
            seed_smiles="CCO",
            n_rounds=1,
        )

    assert any("is the first primary" in r.message for r in caplog.records)


def test_a_tanimoto_annotator_is_scored_against_the_run_s_seed(monkeypatch):
    from chemsampler import execution

    seen = []
    real = execution.LocalExecutor.run_rdkit

    def run_rdkit(self, annotator, smiles_list):
        seen.append(annotator.args())
        return real(self, annotator, smiles_list)

    monkeypatch.setattr(execution.LocalExecutor, "run_rdkit", run_rdkit)
    annotators = [
        AnnotatorSpec("score", StubAnnotator({"CCO": 1.0, "CO": 2.0}), 0.0, "higher"),
        AnnotatorSpec(
            "tanimoto", TanimotoAnnotator(), 0.0, "higher", priority="guard-rail"
        ),
    ]

    _, candidates = hill_climb(
        StubGenerator([["CO"]]), annotators, seed_smiles="CCO", n_rounds=1
    )

    assert seen == [["CCO"], ["CCO"]]  # the seed, then round 1
    assert (
        candidates[1]["tanimoto"].tolist() == candidates[1]["tanimoto_to_seed"].tolist()
    )


# --- annotate ------------------------------------------------------------


def test_annotate_scores_one_smiles_and_reports_cutoffs_satisfied():
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator({"CCO": 1.0}), cutoff=0.0, direction="higher"
        ),
        AnnotatorSpec("missing", StubAnnotator({}), cutoff=0.0, direction="higher"),
    ]

    result = annotate("CCO", annotators)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["source"] == "input"
    assert row["score"] == 1.0
    assert math.isnan(row["missing"])
    assert row["cutoffs_satisfied"] == 1
    assert "tanimoto_to_seed" not in result.columns


def test_annotate_rejects_empty_annotators():
    with pytest.raises(ValueError, match="empty"):
        annotate("CCO", [])


# --- round table sorting -------------------------------------------------


def test_sequential_mode_round_table_sorted_descending_for_higher():
    scores = {"CCO": 1.0, "CCC": 3.0, "CCN": 5.0, "CC": 2.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["score"]) == [5.0, 3.0, 2.0]


def test_sequential_mode_round_table_sorted_ascending_for_lower():
    scores = {"CCO": 1.0, "CCC": 3.0, "CCN": 5.0, "CC": 2.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("inf"), direction="lower"
        )
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["score"]) == [2.0, 3.0, 5.0]


def test_joint_mode_round_table_sorted_by_cutoffs_satisfied_descending():
    a_scores = {"CCO": 1.0, "CCC": 1.0, "CCN": 1.0, "CC": 0.0}
    b_scores = {"CCO": 1.0, "CCC": 0.0, "CCN": 1.0, "CC": 0.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.5, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=0.5, direction="higher"),
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, mode="joint", seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["cutoffs_satisfied"]) == [2, 1, 0]


# --- write_results -----------------------------------------------------


def test_write_results_writes_summary_and_one_csv_per_round(tmp_path):
    summary = pd.DataFrame(
        [
            {
                "round": 0,
                "smiles": "CCO",
                "score": 1.0,
                "is_new_best": True,
                "active_annotator_id": "a",
            }
        ]
    )
    candidates_by_round = {0: pd.DataFrame([{"smiles": "CCO", "source": "seed"}])}

    write_results(summary, candidates_by_round, str(tmp_path / "out"))

    written_summary = pd.read_csv(tmp_path / "out" / "summary.csv")
    assert written_summary["smiles"].tolist() == ["CCO"]
    written_round0 = pd.read_csv(tmp_path / "out" / "round0.csv")
    assert written_round0["source"].tolist() == ["seed"]


def test_write_results_creates_output_dir_including_parents(tmp_path):
    write_results(pd.DataFrame(), {}, str(tmp_path / "a" / "b" / "c"))

    assert (tmp_path / "a" / "b" / "c" / "summary.csv").exists()


# --- similarity -----------------------------------------------------------------


class StubSimilarity:
    """Returns one fixed set of neighbours per call; refuses to run without a seed."""

    model_id = "stub_sim"

    def __init__(self, neighbours):
        self._neighbours = neighbours

    def search_by_model(self, seed_smiles):
        if seed_smiles is None:
            raise SeedRequiredError(f"{self.model_id} requires a seed molecule")
        return {self.model_id: list(self._neighbours)}


def _score_only(scores):
    return [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]


def test_similarity_candidates_join_the_generators_and_keep_provenance():
    annotators = _score_only({"CCO": 1.0, "CCC": 2.0, "CCN": 3.0})

    summary, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        mode="sequential",
        similarity=StubSimilarity(["CCN", "CCC"]),
        seed_smiles="CCO",
        n_rounds=1,
    )

    sources = candidates_by_round[1].set_index("smiles")["source"].to_dict()
    assert sources == {"CCC": "stub|stub_sim", "CCN": "stub_sim"}
    assert summary.iloc[-1]["smiles"] == "CCN"


@pytest.fixture
def verbose_logger():
    """Quiet mode hides warnings raised while generating, and a CLI test can leave it on."""
    previous = logger._quiet_mode
    logger.set_quiet_mode(False)
    yield
    logger.set_quiet_mode(previous)


def test_similarity_is_skipped_with_a_warning_when_there_is_no_seed(
    caplog, verbose_logger
):
    annotators = _score_only({"CCC": 2.0, "CCN": 3.0})

    with caplog.at_level("WARNING"):
        _, candidates_by_round = hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
            mode="sequential",
            similarity=StubSimilarity(["CCN"]),
            n_rounds=1,
        )

    assert candidates_by_round[1]["smiles"].tolist() == ["CCC"]
    assert any(
        "stub_sim" in r.message and "requires a seed" in r.message
        for r in caplog.records
    )


def test_similarity_model_is_rejected_as_generator():
    with pytest.raises(TypeError, match="similarity="):
        hill_climb(
            HubSimilarity("eos8nyg", models_dir="models"),
            _score_only({"CCO": 1.0}),
            mode="sequential",
            seed_smiles="CCO",
        )


# --- guard-rails --------------------------------------------------------------------


def _primary(scores, cutoff=1.5, direction="higher", annotator_id="score"):
    return AnnotatorSpec(
        annotator_id, StubAnnotator(scores), cutoff=cutoff, direction=direction
    )


def _guard(scores, cutoff, direction="higher", annotator_id="guard"):
    return AnnotatorSpec(
        annotator_id,
        StubAnnotator(scores),
        cutoff=cutoff,
        direction=direction,
        priority="guard-rail",
    )


PRIMARY_SCORES = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
GUARD_SCORES = {"CCO": 1.0, "CCC": 1.0, "CCN": 0.0}


@pytest.mark.parametrize(
    "mode, extra",
    [
        ("joint", {}),
        ("weighted", {}),
        ("sequential", {}),
        ("incremental", {"n_steps": 1}),
    ],
)
def test_guard_rail_blocks_the_best_primary_candidate_in_every_mode(mode, extra):
    # CCN is listed first and scores best, so it would win without the guard-rail.
    annotators = [_primary(PRIMARY_SCORES), _guard(GUARD_SCORES, cutoff=0.5)]

    summary, candidates_by_round = hill_climb(
        StubGenerator([["CCN", "CCC"]]),
        annotators,
        mode=mode,
        seed_smiles="CCO",
        n_rounds=1,
        **extra,
    )

    assert summary.iloc[-1]["smiles"] == "CCC"
    round1 = candidates_by_round[1].set_index("smiles")
    assert not round1.loc["CCN", "guardrails_satisfied"]
    assert round1.loc["CCC", "guardrails_satisfied"]


def test_guard_rails_do_not_count_or_weigh():
    annotators = [
        _primary({"CCO": 1.0, "CCC": 2.0}),
        _guard({"CCO": 100.0, "CCC": 100.0}, cutoff=0.0),
    ]

    _, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        mode="weighted",
        seed_smiles="CCO",
        n_rounds=1,
    )

    row = candidates_by_round[1].iloc[0]
    assert row["cutoffs_satisfied"] == 1
    assert row["weighted_score"] == 2.0
    assert row["guard"] == 100.0


def test_guardrails_are_satisfied_when_there_are_none():
    _, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        [_primary(PRIMARY_SCORES)],
        mode="sequential",
        seed_smiles="CCO",
        n_rounds=1,
    )

    assert candidates_by_round[1]["guardrails_satisfied"].all()


def test_range_guard_rail_accepts_only_values_inside_and_never_missing_ones():
    guard = _guard(
        {"CCO": 5.0, "CCC": 5.0, "CCN": 9.0},  # "CCCC" is deliberately unscorable
        cutoff=(4.0, 6.0),
        direction="range",
    )
    annotators = [_primary({"CCO": 1.0, "CCC": 2.0, "CCN": 3.0, "CCCC": 4.0}), guard]

    summary, candidates_by_round = hill_climb(
        StubGenerator([["CCCC", "CCN", "CCC"]]),
        annotators,
        mode="sequential",
        seed_smiles="CCO",
        n_rounds=1,
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert round1["guardrails_satisfied"].to_dict() == {
        "CCCC": False,
        "CCN": False,
        "CCC": True,
    }
    assert summary.iloc[-1]["smiles"] == "CCC"


def test_run_stops_when_every_candidate_breaks_a_guard_rail():
    annotators = [
        _primary(PRIMARY_SCORES),
        _guard({"CCO": 1.0, "CCN": 0.0}, cutoff=0.5),
    ]

    summary, _ = hill_climb(
        StubGenerator([["CCN"]]),
        annotators,
        mode="sequential",
        seed_smiles="CCO",
        n_rounds=3,
    )

    assert list(summary["round"]) == [0]


def test_a_run_needs_at_least_one_primary_annotator():
    with pytest.raises(ValueError, match="primary"):
        hill_climb(
            StubGenerator([["CCC"]]),
            [_guard(GUARD_SCORES, cutoff=0.5)],
            mode="joint",
            seed_smiles="CCO",
        )


def test_guard_rails_never_become_the_active_annotator():
    annotators = [_guard(GUARD_SCORES, cutoff=0.5), _primary(PRIMARY_SCORES)]

    summary, _ = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        mode="sequential",
        seed_smiles="CCO",
        n_rounds=1,
    )

    assert set(summary["active_annotator_id"]) == {"score"}


def test_a_seed_that_breaks_a_guard_rail_only_warns(caplog, verbose_logger):
    annotators = [
        _primary(PRIMARY_SCORES),
        _guard({"CCO": 0.0, "CCC": 1.0}, cutoff=0.5),
    ]

    with caplog.at_level("WARNING"):
        summary, _ = hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
            mode="sequential",
            seed_smiles="CCO",
            n_rounds=1,
        )

    assert any(
        "guard-rail" in r.message and "guard" in r.message for r in caplog.records
    )
    assert summary.iloc[-1]["smiles"] == "CCC"


def test_annotate_accepts_guard_rails_alone_and_reports_them():
    result = annotate("CCO", [_guard({"CCO": 1.0}, cutoff=0.5)])

    row = result.iloc[0]
    assert row["cutoffs_satisfied"] == 0
    assert row["guardrails_satisfied"]


# --- a similarity model without a generator -----------------------------------------


def test_a_similarity_model_can_drive_a_run_alone():
    annotators = _score_only({"CCO": 1.0, "CCC": 2.0, "CCN": 3.0})

    summary, candidates_by_round = hill_climb(
        None,
        annotators,
        mode="sequential",
        similarity=StubSimilarity(["CCN", "CCC"]),
        seed_smiles="CCO",
        n_rounds=1,
    )

    assert set(candidates_by_round[1]["source"]) == {"stub_sim"}
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_a_run_needs_a_generator_or_a_similarity_model():
    with pytest.raises(ValueError, match="generator, a similarity model"):
        hill_climb(None, _score_only({"CCO": 1.0}), mode="joint", seed_smiles="CCO")


# --- a run stopped midway keeps its finished rounds -------------------------------


class StoppingGenerator(StubGenerator):
    """Like StubGenerator, but raises `error` once its rounds of candidates run out."""

    def __init__(self, candidates_by_round, error):
        super().__init__(candidates_by_round)
        self._error = error

    def generate_by_model(self, seed_smiles):
        if self._call_count == len(self._candidates_by_round):
            raise self._error
        return super().generate_by_model(seed_smiles)


def _higher(annotator_id, scores, cutoff=float("-inf")):
    return AnnotatorSpec(
        annotator_id, StubAnnotator(scores), cutoff=cutoff, direction="higher"
    )


@pytest.mark.parametrize(
    "error", [KeyboardInterrupt(), RuntimeError("stub: run.sh failed")]
)
def test_a_run_stopped_midway_keeps_the_rounds_it_finished(error):
    generator = StoppingGenerator([["CCC"]], error)

    with pytest.raises(type(error)) as raised:
        hill_climb(
            generator,
            [_higher("score", {"CCO": 1.0, "CCC": 2.0})],
            mode="weighted",
            seed_smiles="CCO",
            n_rounds=5,
        )

    summary, candidates_by_round = raised.value.partial
    assert list(summary["round"]) == [0, 1]
    assert summary.iloc[-1]["smiles"] == "CCC"
    assert set(candidates_by_round) == {0, 1}


def test_a_stop_in_a_later_stage_keeps_the_earlier_stages():
    # stage a meets its cutoff in round 1; stage b stops in round 2
    annotators = [
        _higher("a", {"CCO": 0.0, "CCC": 2.0}, cutoff=1.0),
        _higher("b", {"CCO": 0.0, "CCC": 0.0}, cutoff=1.0),
    ]

    with pytest.raises(KeyboardInterrupt) as raised:
        hill_climb(
            StoppingGenerator([["CCC"]], KeyboardInterrupt()),
            annotators,
            mode="sequential",
            seed_smiles="CCO",
        )

    summary, candidates_by_round = raised.value.partial
    assert list(summary["round"]) == [0, 1]
    assert list(summary["active_annotator_id"]) == ["a", "a"]
    assert set(candidates_by_round) == {0, 1}


def test_a_stop_within_an_incremental_step_keeps_its_rounds_and_their_level():
    # targets 1.0 then 2.0: round 1 clears step 1, round 2 falls short in step 2
    scores = {"CCO": 0.0, "CCC": 1.0, "CCN": 1.5}

    with pytest.raises(KeyboardInterrupt) as raised:
        hill_climb(
            StoppingGenerator([["CCC"], ["CCN"]], KeyboardInterrupt()),
            [_higher("score", scores, cutoff=2.0)],
            mode="incremental",
            seed_smiles="CCO",
            n_steps=2,
            n_rounds=5,
        )

    summary, _ = raised.value.partial
    assert list(summary["round"]) == [0, 1, 2]
    assert list(summary["step_level"]) == [0, 1, 2]


def test_a_stop_while_scoring_the_seed_has_no_rounds_to_keep():
    class InterruptedAnnotator:
        def score(self, smiles_list):
            raise KeyboardInterrupt

    annotators = [
        AnnotatorSpec("score", InterruptedAnnotator(), cutoff=0.0, direction="higher")
    ]

    with pytest.raises(KeyboardInterrupt) as raised:
        hill_climb(
            StubGenerator([["CCC"]]), annotators, mode="joint", seed_smiles="CCO"
        )

    assert not hasattr(raised.value, "partial")


# --- a stopped run continues from the round tables it saved -----------------------

#: Every round brings one better molecule (awkward floats, to test exact reads).
RESUME_A = {"CCO": 0.1, "CCC": 1 + 1 / 3, "CCN": 2.5, "CCCC": 3.7, "CCCN": 4.9}
RESUME_B = {"CCO": 0.1 + 0.2, "CCC": 1 / 7, "CCN": 3.5, "CCCC": 4.0, "CCCN": 4.1}
RESUME_ROUNDS = [["CCC"], ["CCN"], ["CCCC"], ["CCCN"]]

#: Per mode: cutoffs of a and b (None: b unused) and the round to stop in. Each
#: stop is one the next round depends on: sequential stops in stage b's first
#: round, incremental in step 2's.
RESUME_CASES = {
    "weighted": {"cutoffs": (float("-inf"), float("-inf")), "stop": 3},
    "joint": {"cutoffs": (1.0, 3.0), "stop": 3},
    "sequential": {"cutoffs": (2.0, 3.0), "stop": 3},
    "incremental": {"cutoffs": (3.0, None), "stop": 3, "n_steps": 2},
}


class CountingAnnotator(StubAnnotator):
    def __init__(self, scores):
        super().__init__(scores)
        self.calls = []

    def score(self, smiles_list):
        self.calls.append(list(smiles_list))
        return super().score(smiles_list)


def _resume_annotators(case):
    cutoff_a, cutoff_b = RESUME_CASES[case]["cutoffs"]
    specs = [
        AnnotatorSpec("a", CountingAnnotator(RESUME_A), cutoff_a, direction="higher")
    ]
    if cutoff_b is not None:
        specs.append(
            AnnotatorSpec(
                "b", CountingAnnotator(RESUME_B), cutoff_b, direction="higher"
            )
        )
    return specs


def _climb(case, generator, annotators, checkpoint):
    return hill_climb(
        generator,
        annotators,
        mode=case,
        seed_smiles="CCO",
        n_rounds=4,
        n_steps=RESUME_CASES[case].get("n_steps"),
        checkpoint=checkpoint,
    )


@pytest.mark.parametrize("case", list(RESUME_CASES))
def test_a_resumed_run_ends_exactly_as_an_uninterrupted_one(case, tmp_path):
    stop = RESUME_CASES[case]["stop"]
    whole, whole_rounds = _climb(
        case, StubGenerator(RESUME_ROUNDS), _resume_annotators(case), None
    )
    with pytest.raises(KeyboardInterrupt):
        stopping = StoppingGenerator(RESUME_ROUNDS[: stop - 1], KeyboardInterrupt())
        _climb(case, stopping, _resume_annotators(case), tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        f"round{n}.csv" for n in range(stop)
    ]

    generator = StubGenerator(RESUME_ROUNDS[stop - 1 :])
    annotators = _resume_annotators(case)
    summary, rounds = _climb(case, generator, annotators, tmp_path)

    pd.testing.assert_frame_equal(summary, whole, check_exact=True)
    assert list(rounds) == list(whole_rounds)
    for n, table in rounds.items():
        pd.testing.assert_frame_equal(table, whole_rounds[n], check_exact=True)
    # the rounds it had were read back: only the later ones ran models
    live = len(whole_rounds) - stop
    assert live >= 1 and generator._call_count == live
    assert len(annotators[0].annotator.calls) == live


def test_a_finished_run_called_again_runs_no_model(tmp_path):
    first = _climb(
        "weighted",
        StubGenerator(RESUME_ROUNDS),
        _resume_annotators("weighted"),
        tmp_path,
    )

    again = _climb(
        "weighted", StubGenerator([]), _resume_annotators("weighted"), tmp_path
    )

    pd.testing.assert_frame_equal(again[0], first[0], check_exact=True)


def test_a_half_written_round_is_run_again(tmp_path):
    with pytest.raises(KeyboardInterrupt):
        _climb(
            "weighted",
            StoppingGenerator(RESUME_ROUNDS[:1], KeyboardInterrupt()),
            _resume_annotators("weighted"),
            tmp_path,
        )
    (tmp_path / "round2.csv.partial").write_text(
        "smiles,source\nCC"
    )  # killed mid-write

    generator = StubGenerator(RESUME_ROUNDS[1:])
    summary, _ = _climb("weighted", generator, _resume_annotators("weighted"), tmp_path)

    assert generator._call_count == 3  # rounds 2-4, round 2 included
    assert list(summary["round"]) == [0, 1, 2, 3, 4]


def test_a_saved_round_from_other_annotators_is_refused(tmp_path):
    _climb(
        "weighted",
        StubGenerator(RESUME_ROUNDS),
        _resume_annotators("weighted"),
        tmp_path,
    )

    with pytest.raises(ValueError, match="round0.csv was not made by this run"):
        _climb(
            "weighted", StubGenerator([]), _resume_annotators("incremental"), tmp_path
        )


def test_a_round_killed_while_being_written_leaves_no_table(tmp_path, monkeypatch):
    real = pd.DataFrame.to_csv

    def killed_writing_round_2(self, path, *args, **kwargs):
        if "round2" in str(path):
            with open(path, "w") as f:
                f.write("smiles,source\nCC")
            raise KeyboardInterrupt
        return real(self, path, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", killed_writing_round_2)
    with pytest.raises(KeyboardInterrupt):
        _climb(
            "weighted",
            StubGenerator(RESUME_ROUNDS),
            _resume_annotators("weighted"),
            tmp_path,
        )

    assert not (tmp_path / "round2.csv").exists()
