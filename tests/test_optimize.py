import math

import pandas as pd
import pytest

from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import annotate, hill_climb, write_results


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


def test_tanimoto_cutoff_requires_seed_raises():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "score",
            StubAnnotator({"CCC": 1.0}),
            cutoff=float("-inf"),
            direction="higher",
        )
    ]

    with pytest.raises(ValueError, match="tanimoto_cutoff"):
        hill_climb(
            generator, annotators, mode="joint", seed_smiles=None, tanimoto_cutoff=0.5
        )


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


def test_hill_climb_requires_mode():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "score",
            StubAnnotator({"CCC": 1.0}),
            cutoff=float("-inf"),
            direction="higher",
        )
    ]

    with pytest.raises(TypeError):
        hill_climb(generator, annotators, seed_smiles="CCO")


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


def test_sequential_mode_floor_is_the_achieved_value_not_the_original_cutoff():
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

    # "a" achieves 10.0 in stage 1 (well above its cutoff of 0.0), which becomes
    # the floor for stage 2. CCN clears the *original* cutoff (2.0 >= 0.0) but
    # not the achieved floor, so it can't win despite having the higher "b"
    # score (8.0 > 6.0) - CC must win instead.
    assert summary.iloc[-1]["smiles"] == "CC"


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

    assert list(summary["round"]) == [0]
    assert candidates_by_round[1].iloc[0]["cutoffs_satisfied"] == 0


def test_sequential_stage_with_no_eligible_candidate_stops_entire_run():
    a_scores = {"CCO": 1.0, "CCC": 2.0}
    b_scores = {"CCO": 1.0, "CCC": 2.0}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=0.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator, annotators, mode="sequential", seed_smiles="CCO", n_rounds=1
    )

    # Stage "a" never clears its own cutoff, so stage "b" never starts.
    assert list(summary["active_annotator_id"]) == ["a"]


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
    # cutoff (10.0 >= 0.0), so that value locks in as the floor and stage "b"
    # still runs. (round 0 = seed, round 1 = "a"'s failed attempt, round 2 = "b".)
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


def test_tanimoto_cutoff_excluded_from_joint_mode_count():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    a_scores = {"CCO": 1.0, "CO": 1.0, "c1ccc2ccccc2c1": 1.0}
    b_scores = {"CCO": 1.0, "CO": 0.0, "c1ccc2ccccc2c1": 10.0}
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=0.0, direction="higher"),
        AnnotatorSpec("b", StubAnnotator(b_scores), cutoff=5.0, direction="higher"),
    ]

    summary, candidates_by_round = hill_climb(
        generator,
        annotators,
        mode="joint",
        seed_smiles="CCO",
        n_rounds=1,
        tanimoto_cutoff=0.2,
    )

    # Naphthalene clears both cutoffs (count=2) but fails the Tanimoto gate
    # (0.0 < 0.2); methanol clears only "a" (count=1) but passes Tanimoto
    # (0.2857 >= 0.2) and wins despite the lower count.
    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["c1ccc2ccccc2c1", "cutoffs_satisfied"] == 2
    assert summary.iloc[-1]["smiles"] == "CO"
    assert summary.iloc[-1]["score"] == 1


def test_tanimoto_direction_lower_seeks_novelty():
    # Real Tanimoto (Morgan, radius 2) to "CCO": "CO" = 0.2857, naphthalene = 0.0.
    scores = {"CCO": 1.0, "CO": 1.0, "c1ccc2ccccc2c1": 1.0}
    generator = StubGenerator([["CO", "c1ccc2ccccc2c1"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(scores), cutoff=0.0, direction="higher"),
    ]

    summary, _ = hill_climb(
        generator,
        annotators,
        mode="joint",
        seed_smiles="CCO",
        n_rounds=1,
        tanimoto_cutoff=0.2,
        tanimoto_direction="lower",
    )

    # Both candidates clear the annotator cutoff equally, so the Tanimoto gate
    # decides the winner. Methanol (0.2857) is too similar under a "lower"
    # (novelty) gate and is excluded; naphthalene (0.0) is admitted and wins -
    # the opposite outcome from the same cutoff under the default "higher" gate
    # (see test_tanimoto_cutoff_excluded_from_joint_mode_count above).
    assert summary.iloc[-1]["smiles"] == "c1ccc2ccccc2c1"


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
