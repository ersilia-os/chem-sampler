import pandas as pd
import pytest

from chemsampler.models.hub_annotator import HubAnnotator


class FakeHubModel:
    """Stands in for HubModel, recording each run() call's chunk size."""

    def __init__(self, model_id):
        self.model_id = model_id
        self.calls = []

    def run(self, smiles_list):
        self.calls.append(list(smiles_list))
        return pd.DataFrame(
            {"smiles": smiles_list, "score": [float(len(s)) for s in smiles_list]}
        )


@pytest.fixture
def annotator(monkeypatch):
    monkeypatch.setattr("chemsampler.models.hub_annotator.HubModel", FakeHubModel)
    return HubAnnotator("fake-model")


def test_scores_every_molecule(annotator):
    scores = annotator.score(["C", "CC", "CCC"])

    assert scores == {"C": 1.0, "CC": 2.0, "CCC": 3.0}
    assert len(annotator._hub_model.calls) == 1


def test_chunks_large_inputs(annotator, monkeypatch):
    monkeypatch.setattr("chemsampler.models.hub_annotator._MAX_CHUNK", 3)

    scores = annotator.score(["A", "B", "C", "D", "E", "F", "G"])

    assert len(scores) == 7
    assert [len(c) for c in annotator._hub_model.calls] == [3, 3, 1]


def test_empty_input_makes_no_call(annotator):
    scores = annotator.score([])

    assert scores == {}
    assert annotator._hub_model.calls == []


def test_null_scores_are_dropped_and_warned(annotator, caplog):
    class NullableFakeHubModel(FakeHubModel):
        def run(self, smiles_list):
            df = super().run(smiles_list)
            df.loc[df["smiles"] == "bad", "score"] = None
            return df

    annotator._hub_model = NullableFakeHubModel("fake-model")

    with caplog.at_level("WARNING"):
        scores = annotator.score(["ok", "bad"])

    assert scores == {"ok": 2.0}
    assert any("1 of 2" in r.message for r in caplog.records)
