import csv

import pytest

from chemsampler.config import load_annotators, load_generators
from chemsampler.models.chembl import ChemblSampler
from chemsampler.models.generator import HubGenerator


class FakeHubModel:
    """Stands in for HubModel so constructing a HubGenerator/HubAnnotator is free."""

    def __init__(self, model_id, backend="ersilia"):
        self.model_id = model_id
        self.backend = backend


@pytest.fixture(autouse=True)
def no_real_hub_calls(monkeypatch):
    monkeypatch.setattr("chemsampler.models.generator.HubModel", FakeHubModel)
    monkeypatch.setattr("chemsampler.models.hub_annotator.HubModel", FakeHubModel)


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def test_load_generators_default_excludes_chembl():
    pool = load_generators()

    assert all(not isinstance(g, ChemblSampler) for g in pool.generators)
    assert all(isinstance(g, HubGenerator) for g in pool.generators)


def test_load_generators_includes_chembl_when_explicitly_listed(tmp_path):
    path = tmp_path / "generators.csv"
    _write_csv(path, ["generator_id"], [["eos9taz"], ["chembl"]])

    pool = load_generators(str(path))

    assert {type(g).__name__ for g in pool.generators} == {
        "HubGenerator",
        "ChemblSampler",
    }


def test_load_generators_rejects_missing_column(tmp_path):
    path = tmp_path / "generators.csv"
    _write_csv(path, ["wrong_column"], [["eos9taz"]])

    with pytest.raises(ValueError, match="generator_id"):
        load_generators(str(path))


def test_load_annotators_builds_specs_from_csv(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(
        path,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.3", "higher"], ["eos4zfy", "50", "lower"]],
    )

    specs = load_annotators(str(path))

    assert [spec.annotator_id for spec in specs] == ["qed", "eos4zfy"]
    assert specs[0].cutoff == 0.3
    assert specs[0].direction == "higher"
    assert specs[1].cutoff == 50.0
    assert specs[1].direction == "lower"


def test_load_annotators_preserves_csv_row_order(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(
        path,
        ["annotator_id", "cutoff", "direction"],
        [
            ["eos4zfy", "0", "higher"],
            ["qed", "0.3", "higher"],
            ["eos6ost", "1", "higher"],
        ],
    )

    specs = load_annotators(str(path))

    assert [spec.annotator_id for spec in specs] == ["eos4zfy", "qed", "eos6ost"]


def test_load_annotators_rejects_missing_cutoff(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(path, ["annotator_id", "cutoff", "direction"], [["qed", "", "higher"]])

    with pytest.raises(ValueError, match="cutoff"):
        load_annotators(str(path))


def test_load_annotators_rejects_invalid_cutoff_value(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(
        path,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "not-a-number", "higher"]],
    )

    with pytest.raises(ValueError, match="cutoff"):
        load_annotators(str(path))


def test_load_annotators_rejects_unknown_direction(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(
        path, ["annotator_id", "cutoff", "direction"], [["qed", "0.3", "sideways"]]
    )

    with pytest.raises(ValueError, match="direction"):
        load_annotators(str(path))


def test_load_annotators_rejects_duplicate_annotator_id(tmp_path):
    path = tmp_path / "annotators.csv"
    _write_csv(
        path,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.3", "higher"], ["qed", "0.5", "higher"]],
    )

    with pytest.raises(ValueError, match="duplicate"):
        load_annotators(str(path))
