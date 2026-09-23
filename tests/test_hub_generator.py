import pandas as pd
import pytest

from chemsampler.models.generator import HubGenerator, SeedRequiredError


class FakeHubModel:
    """Stands in for HubModel, returning two fixed candidates per call."""

    def __init__(self, model_id, backend="ersilia"):
        self.model_id = model_id
        self.backend = backend

    def run(self, smiles_list):
        return pd.DataFrame({"smi_0": ["CCN"], "smi_1": ["CCO"]})


@pytest.fixture
def generator(monkeypatch):
    monkeypatch.setattr("chemsampler.models.generator.HubModel", FakeHubModel)
    return HubGenerator("fake-model")


def test_raises_when_seed_is_none(generator):
    with pytest.raises(SeedRequiredError):
        generator.generate(None)


def test_generate_by_model_wraps_generate(generator):
    by_model = generator.generate_by_model("CC")

    assert list(by_model) == ["fake-model"]
    assert sorted(by_model["fake-model"]) == ["CCN", "CCO"]
