from chemsampler.optimize import hill_climb


class StubGenerator:
    """Returns one fixed candidate per seed, regardless of input."""

    def __init__(self, candidates_by_round):
        self._candidates_by_round = list(candidates_by_round)
        self._call_count = 0

    def generate(self, seed_smiles):
        candidates = self._candidates_by_round[self._call_count]
        self._call_count += 1
        return candidates


class StubAnnotator:
    """Scores SMILES by a fixed lookup table."""

    def __init__(self, scores):
        self._scores = scores

    def score(self, smiles_list):
        return {smi: self._scores[smi] for smi in smiles_list}


def test_stops_when_a_round_does_not_improve():
    scores = {"seed": 1.0, "better": 2.0, "worse": 0.5}
    generator = StubGenerator([["better"], ["worse"]])
    annotator = StubAnnotator(scores)

    df = hill_climb("seed", generator, annotator, n_rounds=5)

    assert list(df["round"]) == [0, 1, 2]
    assert list(df["is_new_best"]) == [True, True, False]
    assert df.iloc[-1]["smiles"] == "worse"


def test_respects_n_rounds_cap():
    scores = {"seed": 1.0, "r1": 2.0, "r2": 3.0}
    generator = StubGenerator([["r1"], ["r2"]])
    annotator = StubAnnotator(scores)

    df = hill_climb("seed", generator, annotator, n_rounds=2)

    assert list(df["round"]) == [0, 1, 2]
    assert df["is_new_best"].all()


def test_returns_expected_columns():
    generator = StubGenerator([["seed"]])
    annotator = StubAnnotator({"seed": 1.0})

    df = hill_climb("seed", generator, annotator, n_rounds=1)

    assert list(df.columns) == ["round", "smiles", "score", "is_new_best"]


def test_best_score_trajectory_is_non_decreasing():
    scores = {"seed": 1.0, "r1": 1.5, "r2": 1.2, "r3": 2.0}
    generator = StubGenerator([["r1"], ["r2"], ["r3"]])
    annotator = StubAnnotator(scores)

    df = hill_climb("seed", generator, annotator, n_rounds=5)

    best_so_far = df[df["is_new_best"]]["score"]
    assert list(best_so_far) == sorted(best_so_far)
