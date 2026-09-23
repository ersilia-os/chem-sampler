import csv
import gzip

import pytest

from chemsampler.models.chembl import ChemblSampler


@pytest.fixture
def reference_set(tmp_path):
    """A miniature stand-in for the ChEMBL reference set."""
    path = tmp_path / "chembl_test.csv.gz"
    with gzip.open(path, "wt", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chembl_id", "smiles", "mw"])
        for i in range(50):
            writer.writerow([f"CHEMBL{i}", "C" * (i % 8 + 1), "300.00"])
    return str(path)


def test_draws_requested_number(reference_set):
    sampler = ChemblSampler(n=10, random_state=0, path=reference_set)

    assert len(sampler.generate("CCO")) == 10


def test_random_state_is_reproducible(reference_set):
    a = ChemblSampler(n=10, random_state=42, path=reference_set).generate("CCO")
    b = ChemblSampler(n=10, random_state=42, path=reference_set).generate("CCO")

    assert a == b


def test_seed_is_ignored(reference_set):
    sampler = ChemblSampler(n=10, random_state=42, path=reference_set)

    assert sampler.generate("CCO") == sampler.generate("c1ccccc1")


def test_rejects_sample_larger_than_reference_set(reference_set):
    sampler = ChemblSampler(n=999, path=reference_set)

    with pytest.raises(ValueError, match="reference set holds"):
        sampler.generate("CCO")


def test_generate_by_model_wraps_generate(reference_set):
    sampler = ChemblSampler(n=10, random_state=0, path=reference_set)

    by_model = sampler.generate_by_model("CCO")

    assert list(by_model) == ["chembl"]
    assert by_model["chembl"] == sampler.generate("CCO")


def test_generate_accepts_none_seed(reference_set):
    sampler = ChemblSampler(n=10, random_state=0, path=reference_set)

    assert len(sampler.generate(None)) == 10


def test_pools_with_hub_generators(reference_set):
    from chemsampler.models.generator import GeneratorPool

    class StubGenerator:
        model_id = "stub"

        def generate(self, seed_smiles):
            return ["CCO"]

    pool = GeneratorPool(
        [StubGenerator(), ChemblSampler(n=5, random_state=0, path=reference_set)]
    )
    by_model = pool.generate_by_model("CCO")

    assert set(by_model) == {"stub", "chembl"}
    assert len(by_model["chembl"]) == 5
