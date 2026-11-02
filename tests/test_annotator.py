import re

import pandas as pd
import pytest

from chemsampler.models.annotator import MolecularWeightAnnotator, QEDAnnotator


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


@pytest.mark.parametrize(
    "annotator, name", [(MolecularWeightAnnotator(), "mw"), (QEDAnnotator(), "qed")]
)
def test_an_in_house_annotator_logs_how_many_it_scored_and_how_long(
    annotator, name, caplog
):
    with caplog.at_level("INFO"):
        annotator.score(["CCO", "not_a_smiles"])

    scored = rf"{name}: scored 2 molecules \(\d+\.\ds\)"
    assert any(re.fullmatch(scored, r.message) for r in caplog.records)


def test_the_job_script_gives_the_same_values_as_a_local_score(tmp_path):
    from chemsampler.models import _rdkit_scores

    smiles = ["CCO", "not_a_smiles", "c1ccccc1"]
    pd.DataFrame({"smiles": smiles}).to_csv(tmp_path / "input.csv", index=False)

    _rdkit_scores.main(["qed", str(tmp_path / "input.csv"), str(tmp_path / "out.csv")])

    assert pd.read_csv(tmp_path / "out.csv")["smiles"].tolist() == smiles
    values = _rdkit_scores.read(tmp_path / "out.csv")
    assert values[1] is None
    # exactly equal: the values survive the round trip through the file
    local = QEDAnnotator().score(smiles)
    assert {s: v for s, v in zip(smiles, values) if v is not None} == local
