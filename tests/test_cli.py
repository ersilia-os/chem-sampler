import csv

import pandas as pd
import pytest
from click.testing import CliRunner

from chemsampler.cli import cli


class FakeHubModel:
    """Stands in for HubModel so a fake Hub generator needs no real Ersilia calls."""

    def __init__(self, model_id, backend="ersilia"):
        self.model_id = model_id
        self.backend = backend

    def run(self, smiles_list):
        return pd.DataFrame({"smi_0": ["CCN"]})


@pytest.fixture(autouse=True)
def no_real_hub_calls(monkeypatch):
    monkeypatch.setattr("chemsampler.models.generator.HubModel", FakeHubModel)
    monkeypatch.setattr("chemsampler.models.hub_annotator.HubModel", FakeHubModel)


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def test_run_writes_results_and_reports_best_candidate(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(
        annotators_csv,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.0", "higher"]],
    )
    generators_csv = tmp_path / "generators.csv"
    _write_csv(generators_csv, ["generator_id"], [["eos9taz"]])
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--annotators",
            str(annotators_csv),
            "--generators",
            str(generators_csv),
            "--mode",
            "sequential",
            "--seed-smiles",
            "CCO",
            "--n-rounds",
            "1",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Best candidate" in result.output
    assert (output_dir / "summary.csv").exists()
    assert (output_dir / "round1.csv").exists()


def test_run_requires_mode(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(
        annotators_csv,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.0", "higher"]],
    )

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--annotators",
            str(annotators_csv),
            "--seed-smiles",
            "CCO",
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "--mode" in result.output


def test_run_reports_value_error_in_red_and_exits_nonzero(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(annotators_csv, ["annotator_id", "cutoff", "direction"], [])
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--annotators",
            str(annotators_csv),
            "--mode",
            "joint",
            "--seed-smiles",
            "CCO",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 1
    assert not output_dir.exists()


def test_run_handles_empty_summary_without_seed(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(
        annotators_csv,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.0", "higher"]],
    )
    generators_csv = tmp_path / "generators.csv"
    _write_csv(generators_csv, ["generator_id"], [["eos9taz"]])
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--annotators",
            str(annotators_csv),
            "--generators",
            str(generators_csv),
            "--mode",
            "sequential",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "No candidates" in result.output
    assert (output_dir / "summary.csv").exists()


def test_annotate_prints_values_and_cutoff_summary(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(
        annotators_csv,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.0", "higher"]],
    )

    result = CliRunner().invoke(
        cli,
        ["annotate", "--annotators", str(annotators_csv), "--smiles", "CCO"],
    )

    assert result.exit_code == 0, result.output
    assert "qed" in result.output
    assert "cutoffs satisfied" in result.output


def test_annotate_requires_smiles(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(
        annotators_csv,
        ["annotator_id", "cutoff", "direction"],
        [["qed", "0.0", "higher"]],
    )

    result = CliRunner().invoke(cli, ["annotate", "--annotators", str(annotators_csv)])

    assert result.exit_code != 0
    assert "--smiles" in result.output


def test_annotate_reports_value_error_in_red_and_exits_nonzero(tmp_path):
    annotators_csv = tmp_path / "annotators.csv"
    _write_csv(annotators_csv, ["annotator_id", "cutoff", "direction"], [])

    result = CliRunner().invoke(
        cli,
        ["annotate", "--annotators", str(annotators_csv), "--smiles", "CCO"],
    )

    assert result.exit_code == 1
