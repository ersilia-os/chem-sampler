import math

import pandas as pd
import pytest

from chemsampler.models.annotator import TanimotoAnnotator
from chemsampler.models.generator import SeedRequiredError
from chemsampler.models.similarity import HubSimilarity
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import annotate, hill_climb, write_results
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


# --- rounds and stopping --------------------------------------------------------


def test_stops_when_a_round_does_not_improve():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 0.5}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

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
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=2)

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

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=1)

    assert list(summary.columns) == [
        "round",
        "smiles",
        "score",
        "is_new_best",
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
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1].columns) == [
        "smiles",
        "source",
        "score",
        "cutoff_satisfied",
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

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=5)

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

    summary, candidates_by_round = hill_climb(generator, annotators, n_rounds=1)

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


def test_a_run_maximizes_its_primary():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=2)

    assert list(summary["is_new_best"]) == [True, True, True]
    assert summary.iloc[-1]["smiles"] == "CCN"
    assert summary.iloc[-1]["score"] == 3.0


# --- validation ------------------------------------------------------------


def test_hill_climb_rejects_empty_annotators():
    with pytest.raises(ValueError, match="empty"):
        hill_climb(StubGenerator([]), [], seed_smiles="CCO")


def test_hill_climb_rejects_reserved_annotator_id():
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec(
            "cutoff_satisfied",
            StubAnnotator({"CCO": 1.0}),
            cutoff=0.0,
            direction="higher",
        )
    ]

    with pytest.raises(ValueError, match="reserved"):
        hill_climb(generator, annotators, seed_smiles="CCO")


def test_a_run_takes_exactly_one_primary():
    scores = {"CCO": 0.0, "CCC": 1.0}
    two = [
        AnnotatorSpec("a", StubAnnotator(scores), 0.0, "higher"),
        AnnotatorSpec("b", StubAnnotator(scores), 0.0, "higher"),
    ]

    # refused before anything is generated: StubGenerator([]) fails if called
    with pytest.raises(ValueError, match=r"one primary annotator only, got 2 \(a, b\)"):
        hill_climb(StubGenerator([]), two, seed_smiles="CCO")
    with pytest.raises(ValueError, match="one primary annotator only"):
        annotate("CCO", two)


# --- the primary's cutoff --------------------------------------------------------


def test_an_improvement_counts_before_the_cutoff_is_met():
    a_scores = {"CCO": 1.0, "CCC": 2.0}  # CCC beats the seed but not the cutoff
    generator = StubGenerator([["CCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(a_scores), cutoff=5.0, direction="higher")
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    # the improvement (1.0 -> 2.0) is recorded although the cutoff (5.0) is not met
    assert list(summary["round"]) == [0, 1]
    assert summary.iloc[1]["score"] == 2.0
    assert candidates_by_round[1].iloc[0]["cutoff_satisfied"] == False


def test_meeting_the_cutoff_stops_the_run():
    scores = {"CCO": 1.0, "CCC": 2.0, "CCN": 3.0}
    generator = StubGenerator([["CCC"], ["CCN"]])
    annotators = [AnnotatorSpec("a", StubAnnotator(scores), 2.0, "higher")]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=5
    )

    assert list(summary["round"]) == [0, 1]
    assert candidates_by_round[1]["cutoff_satisfied"].tolist() == [True]


@pytest.mark.parametrize("direction", ["higher", "lower"])
def test_a_primary_without_a_cutoff_is_optimized(direction, caplog, verbose_logger):
    sign = 1 if direction == "higher" else -1
    scores = {"CCO": 0.0, "CCC": 1.0, "CCN": 2.0, "CC": 3.0, "CCCC": 2.5}
    scores = {smi: sign * value for smi, value in scores.items()}
    generator = StubGenerator([["CCC"], ["CCN"], ["CC"], ["CCCC"]])
    annotators = [
        AnnotatorSpec("a", StubAnnotator(scores), float("nan"), direction),
    ]

    with caplog.at_level("WARNING"):
        summary, candidates_by_round = hill_climb(
            generator, annotators, seed_smiles="CCO", n_rounds=10
        )

    # it climbs while a round improves, and stops at the first that does not
    assert list(summary["smiles"]) == ["CCO", "CCC", "CCN", "CC", "CCCC"]
    assert list(summary["is_new_best"]) == [True, True, True, True, False]
    assert all(t["cutoff_satisfied"].isna().all() for t in candidates_by_round.values())
    assert not any("cutoff" in r.message for r in caplog.records)


def test_a_guard_rail_needs_a_cutoff():
    with pytest.raises(ValueError, match="a guard-rail needs a cutoff"):
        AnnotatorSpec(
            "g", StubAnnotator({}), float("nan"), "higher", priority="guard-rail"
        )


def test_a_candidate_the_primary_cannot_score_cannot_win():
    # "CCC" is deliberately absent from "a"'s table.
    a_scores = {"CCO": 1.0, "CCN": 0.5}
    generator = StubGenerator([["CCC", "CCN"]])
    annotators = [
        AnnotatorSpec(
            "a", StubAnnotator(a_scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    summary, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    round1 = candidates_by_round[1].set_index("smiles")
    assert math.isnan(round1.loc["CCC", "a"])
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_a_round_the_primary_cannot_score_at_all_stops_the_run():
    annotators = [AnnotatorSpec("a", StubAnnotator({"CCO": 1.0}), 5.0, "higher")]

    summary, candidates_by_round = hill_climb(
        StubGenerator([["CCC", "CCN"]]), annotators, seed_smiles="CCO", n_rounds=3
    )

    assert list(summary["round"]) == [0]
    assert candidates_by_round[1]["a"].isna().all()


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

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=1)

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
        generator, annotators, seed_smiles="CCO", n_rounds=1
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

    summary, _ = hill_climb(generator, annotators, seed_smiles="CCO", n_rounds=1)

    assert summary.iloc[-1]["smiles"] == "c1ccc2ccccc2c1"
    assert summary.iloc[-1]["score"] == 0.0


def test_a_higher_tanimoto_primary_warns(caplog, verbose_logger):
    annotators = [AnnotatorSpec("tanimoto", TanimotoAnnotator(), 0.5, "higher")]

    with caplog.at_level("WARNING"):
        hill_climb(
            StubGenerator([["CO"], ["CO"]]),
            annotators,
            seed_smiles="CCO",
            n_rounds=1,
        )

    assert any("is the primary" in r.message for r in caplog.records)


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


def test_annotate_scores_one_smiles_and_reports_the_primary_cutoff():
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator({"CCO": 1.0}), cutoff=0.0, direction="higher"
        ),
        AnnotatorSpec(
            "missing", StubAnnotator({}), 0.0, "higher", priority="guard-rail"
        ),
    ]

    result = annotate("CCO", annotators)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["source"] == "input"
    assert row["score"] == 1.0
    assert math.isnan(row["missing"])
    assert row["cutoff_satisfied"] == True
    assert not row["guardrails_satisfied"]  # a missing score never satisfies one
    assert "tanimoto_to_seed" not in result.columns


def test_annotate_rejects_empty_annotators():
    with pytest.raises(ValueError, match="empty"):
        annotate("CCO", [])


# --- round table sorting -------------------------------------------------


def test_round_table_sorted_descending_for_higher():
    scores = {"CCO": 1.0, "CCC": 3.0, "CCN": 5.0, "CC": 2.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("-inf"), direction="higher"
        )
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["score"]) == [5.0, 3.0, 2.0]


def test_round_table_sorted_ascending_for_lower():
    scores = {"CCO": 1.0, "CCC": 3.0, "CCN": 5.0, "CC": 2.0}
    generator = StubGenerator([["CCC", "CCN", "CC"]])
    annotators = [
        AnnotatorSpec(
            "score", StubAnnotator(scores), cutoff=float("inf"), direction="lower"
        )
    ]

    _, candidates_by_round = hill_climb(
        generator, annotators, seed_smiles="CCO", n_rounds=1
    )

    assert list(candidates_by_round[1]["score"]) == [2.0, 3.0, 5.0]


# --- write_results -----------------------------------------------------


def test_write_results_writes_summary_and_one_csv_per_round(tmp_path):
    summary = pd.DataFrame(
        [
            {
                "round": 0,
                "smiles": "CCO",
                "score": 1.0,
                "is_new_best": True,
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

    def __init__(self, neighbours, model_id="stub_sim"):
        self._neighbours = neighbours
        self.model_id = model_id

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


def test_guard_rail_blocks_the_best_primary_candidate():
    # CCN is listed first and scores best, so it would win without the guard-rail.
    annotators = [_primary(PRIMARY_SCORES), _guard(GUARD_SCORES, 0.5)]

    summary, candidates_by_round = hill_climb(
        StubGenerator([["CCN", "CCC"]]),
        annotators,
        seed_smiles="CCO",
        n_rounds=1,
    )

    assert summary.iloc[-1]["smiles"] == "CCC"
    round1 = candidates_by_round[1].set_index("smiles")
    assert not round1.loc["CCN", "guardrails_satisfied"]
    assert round1.loc["CCC", "guardrails_satisfied"]


def test_guard_rails_do_not_count_toward_the_primary_cutoff():
    annotators = [
        _primary({"CCO": 1.0, "CCC": 1.0}),  # below its cutoff, 1.5
        _guard({"CCO": 100.0, "CCC": 100.0}, cutoff=0.0),
    ]

    _, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        annotators,
        seed_smiles="CCO",
        n_rounds=1,
    )

    row = candidates_by_round[1].iloc[0]
    assert row["cutoff_satisfied"] == False  # the primary's, whatever the guard-rail
    assert row["guardrails_satisfied"]
    assert row["guard"] == 100.0


def test_guardrails_are_satisfied_when_there_are_none():
    _, candidates_by_round = hill_climb(
        StubGenerator([["CCC"]]),
        [_primary(PRIMARY_SCORES)],
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
        seed_smiles="CCO",
        n_rounds=3,
    )

    assert list(summary["round"]) == [0]


def test_a_run_needs_at_least_one_primary_annotator():
    with pytest.raises(ValueError, match="primary"):
        hill_climb(
            StubGenerator([["CCC"]]),
            [_guard(GUARD_SCORES, cutoff=0.5)],
            seed_smiles="CCO",
        )


def test_a_seed_that_breaks_a_guard_rail_only_warns(caplog, verbose_logger):
    annotators = [
        _primary(PRIMARY_SCORES),
        _guard({"CCO": 0.0, "CCC": 1.0}, cutoff=0.5),
    ]

    with caplog.at_level("WARNING"):
        summary, _ = hill_climb(
            StubGenerator([["CCC"]]),
            annotators,
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
    assert pd.isna(row["cutoff_satisfied"])  # no primary to judge
    assert row["guardrails_satisfied"]


# --- a similarity model without a generator -----------------------------------------


def test_a_similarity_model_can_drive_a_run_alone():
    annotators = _score_only({"CCO": 1.0, "CCC": 2.0, "CCN": 3.0})

    summary, candidates_by_round = hill_climb(
        None,
        annotators,
        similarity=StubSimilarity(["CCN", "CCC"]),
        seed_smiles="CCO",
        n_rounds=1,
    )

    assert set(candidates_by_round[1]["source"]) == {"stub_sim"}
    assert summary.iloc[-1]["smiles"] == "CCN"


def test_several_similarity_models_all_propose_candidates():
    annotators = _score_only({"CCO": 1.0, "CCC": 2.0, "CCN": 3.0, "CCCl": 4.0})

    summary, candidates_by_round = hill_climb(
        None,
        annotators,
        similarity=[
            StubSimilarity(["CCN"], model_id="sim_a"),
            StubSimilarity(["CCCl", "CCC"], model_id="sim_b"),
        ],
        seed_smiles="CCO",
        n_rounds=1,
    )

    sources = candidates_by_round[1].set_index("smiles")["source"].to_dict()
    assert sources == {"CCN": "sim_a", "CCCl": "sim_b", "CCC": "sim_b"}
    assert summary.iloc[-1]["smiles"] == "CCCl"


def test_a_run_needs_a_generator_or_a_similarity_model():
    with pytest.raises(ValueError, match="generator, a similarity model"):
        hill_climb(None, _score_only({"CCO": 1.0}), seed_smiles="CCO")


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
            seed_smiles="CCO",
            n_rounds=5,
        )

    summary, candidates_by_round = raised.value.partial
    assert list(summary["round"]) == [0, 1]
    assert summary.iloc[-1]["smiles"] == "CCC"
    assert set(candidates_by_round) == {0, 1}


def test_a_stop_while_scoring_the_seed_has_no_rounds_to_keep():
    class InterruptedAnnotator:
        def score(self, smiles_list):
            raise KeyboardInterrupt

    annotators = [
        AnnotatorSpec("score", InterruptedAnnotator(), cutoff=0.0, direction="higher")
    ]

    with pytest.raises(KeyboardInterrupt) as raised:
        hill_climb(StubGenerator([["CCC"]]), annotators, seed_smiles="CCO")

    assert not hasattr(raised.value, "partial")


# --- a stopped run continues from the round tables it saved -----------------------

#: Every round brings one better molecule (awkward floats, to test exact reads).
RESUME_A = {"CCO": 0.1, "CCC": 1 + 1 / 3, "CCN": 2.5, "CCCC": 3.7, "CCCN": 4.9}
RESUME_B = {"CCO": 0.1 + 0.2, "CCC": 1 / 7, "CCN": 3.5, "CCCC": 4.0, "CCCN": 4.1}
RESUME_ROUNDS = [["CCC"], ["CCN"], ["CCCC"], ["CCCN"]]

#: The primary's cutoff (met in round 4, or none) and the round to stop in.
RESUME_CASES = {
    "cutoff": {"cutoff": 4.0, "stop": 3},
    "optimize": {"cutoff": float("nan"), "stop": 3},
}


class CountingAnnotator(StubAnnotator):
    def __init__(self, scores):
        super().__init__(scores)
        self.calls = []

    def score(self, smiles_list):
        self.calls.append(list(smiles_list))
        return super().score(smiles_list)


def _resume_annotators(case):
    cutoff = RESUME_CASES[case]["cutoff"]
    return [
        AnnotatorSpec("a", CountingAnnotator(RESUME_A), cutoff, direction="higher"),
        AnnotatorSpec(
            "b", CountingAnnotator(RESUME_B), 0.0, "higher", priority="guard-rail"
        ),
    ]


def _climb(case, generator, annotators, checkpoint):
    return hill_climb(
        generator,
        annotators,
        seed_smiles="CCO",
        n_rounds=4,
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
        "optimize",
        StubGenerator(RESUME_ROUNDS),
        _resume_annotators("optimize"),
        tmp_path,
    )

    again = _climb(
        "optimize", StubGenerator([]), _resume_annotators("optimize"), tmp_path
    )

    pd.testing.assert_frame_equal(again[0], first[0], check_exact=True)


def test_a_half_written_round_is_run_again(tmp_path):
    with pytest.raises(KeyboardInterrupt):
        _climb(
            "optimize",
            StoppingGenerator(RESUME_ROUNDS[:1], KeyboardInterrupt()),
            _resume_annotators("optimize"),
            tmp_path,
        )
    (tmp_path / "round2.csv.partial").write_text(
        "smiles,source\nCC"
    )  # killed mid-write

    generator = StubGenerator(RESUME_ROUNDS[1:])
    summary, _ = _climb("optimize", generator, _resume_annotators("optimize"), tmp_path)

    assert generator._call_count == 3  # rounds 2-4, round 2 included
    assert list(summary["round"]) == [0, 1, 2, 3, 4]


def test_a_saved_round_from_other_annotators_is_refused(tmp_path):
    _climb(
        "optimize",
        StubGenerator(RESUME_ROUNDS),
        _resume_annotators("optimize"),
        tmp_path,
    )

    other = [_resume_annotators("optimize")[0], _guard(RESUME_B, 0.0, "higher", "c")]
    with pytest.raises(ValueError, match="round0.csv was not made by this run"):
        _climb("optimize", StubGenerator([]), other, tmp_path)


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
            "optimize",
            StubGenerator(RESUME_ROUNDS),
            _resume_annotators("optimize"),
            tmp_path,
        )

    assert not (tmp_path / "round2.csv").exists()
