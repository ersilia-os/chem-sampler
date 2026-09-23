import math

import pytest

from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import hill_climb


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


def test_stops_when_a_round_does_not_improve():
    # With one directing annotator, no controlling annotators and no Tanimoto
    # cutoff, this reproduces the pre-redesign hill_climb's exact behavior.
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 0.5}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [AnnotatorSpec("score", StubAnnotator(scores), role="directing")]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=5
    )

    assert list(summary["round"]) == [0, 1, 2]
    assert list(summary["is_new_best"]) == [True, True, False]
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert set(candidates_by_round) == {0, 1, 2}


def test_respects_n_rounds_cap():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [AnnotatorSpec("score", StubAnnotator(scores), role="directing")]

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=2)

    assert list(summary["round"]) == [0, 1, 2]
    assert summary["is_new_best"].all()


def test_summary_has_expected_columns():
    generator = StubGenerator([["CCO"]])
    annotators = [AnnotatorSpec("score", StubAnnotator({"CCO": 1.0}), role="directing")]

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=1)

    assert list(summary.columns) == ["round", "smiles", "score", "is_new_best"]


def test_candidates_table_has_expected_columns():
    generator = StubGenerator([["CCO"]])
    annotators = [AnnotatorSpec("score", StubAnnotator({"CCO": 1.0}), role="directing")]

    _, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1].columns) == [
        "smiles",
        "source",
        "score",
        "passes_constraints",
        "tanimoto_to_seed",
    ]


def test_best_score_trajectory_is_non_decreasing():
    scores = {"CCO": 1.0, "CCC": 1.5, "CCN": 1.2, "CC": 2.0}
    generator = StubGenerator([["CCC"], ["CCN"], ["CC"]])
    annotators = [AnnotatorSpec("score", StubAnnotator(scores), role="directing")]

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=5)

    best_so_far = summary[summary["is_new_best"]]["score"]
    assert list(best_so_far) == sorted(best_so_far)


def test_candidate_failing_constraint_is_kept_but_cannot_win():
    directing_scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 1.5}
    controlling_scores = {"CCO": 200.0, "CCC": 999.0, "CCN": 100.0}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(directing_scores), role="directing"),
        AnnotatorSpec(
            "mw", StubAnnotator(controlling_scores), role="controlling", min=0, max=500
        ),
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert round1.loc["CCC", "passes_constraints"] == False  # noqa: E712
    assert round1.loc["CCN", "passes_constraints"] == True  # noqa: E712
    # CCC scores higher but fails the constraint, so CCN wins despite the lower score.
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_round_breaks_when_no_candidate_passes_constraints(caplog):
    directing_scores = {"CCO": 1.0, "CCC": 2.0}
    controlling_scores = {"CCO": 200.0, "CCC": 999.0}
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(directing_scores), role="directing"),
        AnnotatorSpec(
            "mw", StubAnnotator(controlling_scores), role="controlling", min=0, max=500
        ),
    ]

    with caplog.at_level("WARNING"):
        summary, candidates_by_round = hill_climb(
            generator, annotators, seed_smiles="CCO", n_rounds=5
        )

    assert list(summary["round"]) == [0]
    assert len(candidates_by_round[1]) == 1
    assert candidates_by_round[1].iloc[0]["passes_constraints"] == False  # noqa: E712
    assert any("no candidate satisfied" in r.message for r in caplog.records)


def test_tolerance_requires_margin_to_count_as_improvement():
    scores = {"CCO": 1.0, "CCC": 1.05}
    generator = StubGenerator([["CCC"]])
    annotators = [AnnotatorSpec("score", StubAnnotator(scores), role="directing")]

    summary, _ = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=5, tolerance=0.1
    )

    assert list(summary["round"]) == [0, 1]
    assert list(summary["is_new_best"]) == [True, False]


def test_seedless_first_round_is_unconditionally_the_new_best():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator({"CCC": -5.0}), role="directing")
    ]

    summary, candidates_by_round = hill_climb(generator, annotators, n_rounds=1)

    assert list(summary["round"]) == [1]
    assert summary.iloc[0]["is_new_best"] == True  # noqa: E712
    assert 0 not in candidates_by_round
    assert "tanimoto_to_seed" not in candidates_by_round[1].columns


def test_tanimoto_cutoff_requires_seed_raises():
    generator = StubGenerator([["CCC"]])
    annotators = [AnnotatorSpec("score", StubAnnotator({"CCC": 1.0}), role="directing")]

    with pytest.raises(ValueError, match="tanimoto_cutoff"):
        hill_climb(generator, annotators, seed_smiles=None, tanimoto_cutoff=0.5)


def test_hill_climb_rejects_annotators_without_exactly_one_directing():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "a",
            StubAnnotator({"CCO": 1.0, "CCC": 1.0}),
            role="controlling",
            min=0,
            max=1,
        ),
    ]

    with pytest.raises(ValueError, match="directing"):
        hill_climb(generator, annotators, seed_smiles="CCO")


def test_candidate_missing_directing_score_cannot_win():
    # "CCC" is deliberately absent from the directing annotator's table.
    directing_scores = {"CCO": 1.0, "CCN": 0.5}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec("score", StubAnnotator(directing_scores), role="directing")
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert math.isnan(round1.loc["CCC", "score"])
    assert summary.iloc[-1]["smiles"] == "CCN"
