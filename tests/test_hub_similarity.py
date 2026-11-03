import re

import pandas as pd
import pytest

from chemsampler.models.generator import SeedRequiredError
from chemsampler.models.similarity import HubSimilarity


class FakeHubModel:
    """Stands in for HubModel, returning ranked neighbours that include the seed."""

    def __init__(self, model_id, models_dir=None, envs_dir=None):
        self.model_id = model_id
        self.elapsed = 0.0

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


def test_a_search_logs_what_it_found_and_how_long_it_took(similarity, caplog):
    with caplog.at_level("INFO"):
        similarity.search("CCO")

    found = r"eos8nyg: found 2 neighbours for this seed \(\d+\.\ds\)"
    assert any(re.fullmatch(found, r.message) for r in caplog.records)


def test_search_by_model_is_keyed_by_model_id(similarity):
    assert similarity.search_by_model("CCO") == {"eos8nyg": ["CCN", "CCC"]}


def test_rejects_a_model_that_is_not_a_similarity_model():
    with pytest.raises(ValueError, match="not meant to be used in chemsampler"):
        HubSimilarity("eos9taz")


@pytest.mark.parametrize(
    "model", ["eos8nyg", "eos1bgt", "eos80al", "eos3pu0", "eos9b1c", "eos8698"]
)
def test_the_allowed_similarity_models_are_accepted(monkeypatch, model):
    monkeypatch.setattr("chemsampler.models.similarity.HubModel", FakeHubModel)

    assert HubSimilarity(model).search_by_model("CCO") == {model: ["CCN", "CCC"]}
