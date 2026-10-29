import pandas as pd
import pytest

from chemsampler.models.generator import SeedRequiredError
from chemsampler.models.similarity import HubSimilarity


class FakeHubModel:
    """Stands in for HubModel, returning ranked neighbours that include the seed."""

    def __init__(self, model_id, backend="ersilia"):
        self.model_id = model_id
        self.backend = backend

    def run(self, smiles_list):
        return pd.DataFrame(
            {
                "smi_00": ["CCO"],
                "smi_01": ["CCN"],
                "smi_02": ["not-a-smiles"],
                "smi_03": ["CCC"],
                "smi_04": ["CCN"],
            }
        )


@pytest.fixture
def similarity(monkeypatch):
    monkeypatch.setattr("chemsampler.models.similarity.HubModel", FakeHubModel)
    return HubSimilarity("eos8nyg")


def test_search_drops_seed_unparsable_and_duplicates_keeping_rank_order(similarity):
    assert similarity.search("CCO") == ["CCN", "CCC"]


def test_raises_when_seed_is_none(similarity):
    with pytest.raises(SeedRequiredError):
        similarity.search(None)


def test_search_by_model_is_keyed_by_model_id(similarity):
    assert similarity.search_by_model("CCO") == {"eos8nyg": ["CCN", "CCC"]}


def test_rejects_a_model_that_is_not_a_similarity_model():
    with pytest.raises(ValueError, match="not a supported similarity model"):
        HubSimilarity("eos9taz")
