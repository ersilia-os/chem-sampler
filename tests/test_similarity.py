import pytest

from chemsampler.utils.similarity import tanimoto_to_seed


def test_identical_molecule_is_tanimoto_1():
    scores = tanimoto_to_seed("CCO", ["CCO"])

    assert scores["CCO"] == pytest.approx(1.0)


def test_dissimilar_is_less_similar_than_identical():
    scores = tanimoto_to_seed("c1ccccc1", ["c1ccccc1", "CCCCCCCC"])

    assert scores["c1ccccc1"] > scores["CCCCCCCC"]


def test_skips_unparseable_candidates():
    scores = tanimoto_to_seed("CCO", ["CCO", "not_a_smiles"])

    assert set(scores) == {"CCO"}


def test_raises_on_unparseable_seed():
    with pytest.raises(ValueError, match="seed"):
        tanimoto_to_seed("not_a_smiles", ["CCO"])
