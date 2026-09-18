from chemsampler.models.annotator import QEDAnnotator


def test_scores_known_smiles():
    annotator = QEDAnnotator()
    benzene = "c1ccccc1"

    scores = annotator.score([benzene])

    assert set(scores) == {benzene}
    assert 0.0 <= scores[benzene] <= 1.0


def test_skips_invalid_smiles():
    annotator = QEDAnnotator()

    scores = annotator.score(["not_a_smiles", "c1ccccc1"])

    assert set(scores) == {"c1ccccc1"}
