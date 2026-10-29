import csv
import json
import re

import pandas as pd
import pytest
from click.testing import CliRunner

from chemsampler import config
from chemsampler.cli import cli
from chemsampler.hub import catalog
from chemsampler.hub.envs import write_marker
from chemsampler.workspace import create_workspace


class FakeHubModel:
    """Stands in for HubModel so fake Hub models need no real Ersilia calls."""

    def __init__(self, model_id, backend="ersilia", models_dir=None, envs_dir=None):
        self.model_id = model_id
        self.backend = backend

    def run(self, smiles_list):
        return pd.DataFrame({"smi_0": ["CCN"]})


@pytest.fixture(autouse=True)
def no_real_hub_calls(monkeypatch):
    for module in ("generator", "hub_annotator", "similarity"):
        monkeypatch.setattr(f"chemsampler.models.{module}.HubModel", FakeHubModel)
    monkeypatch.setattr(config, "validate_models", lambda **kwargs: None)


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


ANNOTATOR_HEADER = ["source", "model", "cutoff", "direction", "priority"]
QED = ["rdkit", "qed", "0.0", "higher", "primary"]


def _inputs(tmp_path, annotators=(QED,)):
    folder = tmp_path / "inputs"
    folder.mkdir()
    _write_csv(folder / "similarity.csv", ["similarity"], [["eos8nyg"]])
    _write_csv(folder / "generators.csv", ["generators"], [["eos9taz"]])
    _write_csv(folder / "annotators.csv", ANNOTATOR_HEADER, list(annotators))
    return folder


def _run_without_output(*args):
    return CliRunner().invoke(
        cli, ["run", "--mode", "sequential", "--seed-smiles", "CCO", *args]
    )


def _run(tmp_path, *args):
    return CliRunner().invoke(
        cli,
        [
            "run",
            "--mode",
            "sequential",
            "--seed-smiles",
            "CCO",
            "--n-rounds",
            "1",
            "--backend",
            "ersilia",
            "--output-dir",
            str(tmp_path / "out"),
            *args,
        ],
    )


# --- run ----------------------------------------------------------------------------


def test_run_writes_results_and_merges_generator_and_similarity_candidates(tmp_path):
    folder = _inputs(tmp_path)

    result = _run(tmp_path, "--path-to-models", str(folder))

    assert result.exit_code == 0, result.output
    assert "Best candidate" in result.output
    assert (tmp_path / "out" / "summary.csv").exists()
    round1 = pd.read_csv(tmp_path / "out" / "round1.csv")
    assert round1["source"].tolist() == ["eos9taz|eos8nyg"]


def test_run_accepts_everything_in_a_config_file(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(
        json.dumps(
            {
                "similarity": ["eos8nyg"],
                "generators": ["eos9taz"],
                "annotators": [
                    {
                        "source": "rdkit",
                        "model": "qed",
                        "cutoff": 0.0,
                        "direction": "higher",
                        "priority": "primary",
                    }
                ],
            }
        )
    )

    result = _run(tmp_path, "--config", str(path))

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "round1.csv").exists()


def test_run_flags_override_the_config_file(tmp_path):
    folder = _inputs(tmp_path)
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"path_to_models": str(tmp_path / "nowhere")}))

    result = _run(tmp_path, "--config", str(path), "--path-to-models", str(folder))

    assert result.exit_code == 0, result.output


def test_a_secondary_annotator_acts_as_a_guard_rail_through_the_cli(tmp_path):
    folder = _inputs(tmp_path, [QED, ["rdkit", "mw", "0-10", "range", "secondary"]])

    result = _run(tmp_path, "--path-to-models", str(folder))

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(tmp_path / "out" / "round1.csv")
    assert round1["guardrails_satisfied"].tolist() == [False]  # CCN weighs about 45 Da


def test_run_requires_mode(tmp_path):
    folder = _inputs(tmp_path)

    result = CliRunner().invoke(
        cli,
        ["run", "--path-to-models", str(folder), "--output-dir", str(tmp_path / "o")],
    )

    assert result.exit_code != 0
    assert "--mode" in result.output


def test_run_needs_the_inputs_and_reports_errors_in_red(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    empty = tmp_path / "empty"
    empty.mkdir()

    result = _run(tmp_path, "--path-to-models", str(empty))

    assert result.exit_code == 1
    assert "give at least one of similarity and generators" in result.output
    assert not (tmp_path / "out").exists()


def test_run_sh_needs_the_repositories_folder(tmp_path):
    folder = _inputs(tmp_path)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--mode",
            "sequential",
            "--seed-smiles",
            "CCO",
            "--path-to-models",
            str(folder),
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 1
    assert "path-to-model-repositories" in result.output


def test_run_handles_empty_summary_without_seed(tmp_path):
    folder = _inputs(tmp_path)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--mode",
            "joint",
            "--path-to-models",
            str(folder),
            "--backend",
            "ersilia",
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "No candidates were generated" in result.output


def test_run_takes_each_input_as_its_own_csv_flag(tmp_path):
    folder = _inputs(tmp_path)

    result = _run(
        tmp_path,
        "--similarity",
        str(folder / "similarity.csv"),
        "--generators",
        str(folder / "generators.csv"),
        "--annotators",
        str(folder / "annotators.csv"),
    )

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(tmp_path / "out" / "round1.csv")
    assert round1["source"].tolist() == ["eos9taz|eos8nyg"]


def test_run_can_use_a_similarity_model_without_generators(tmp_path):
    folder = _inputs(tmp_path)
    (folder / "generators.csv").unlink()

    result = _run(tmp_path, "--path-to-models", str(folder))

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(tmp_path / "out" / "round1.csv")
    assert round1["source"].tolist() == ["eos8nyg"]


def test_run_can_use_generators_without_a_similarity_model(tmp_path):
    folder = _inputs(tmp_path)
    (folder / "similarity.csv").unlink()

    result = _run(tmp_path, "--path-to-models", str(folder))

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(tmp_path / "out" / "round1.csv")
    assert round1["source"].tolist() == ["eos9taz"]


def test_run_finds_its_inputs_in_the_working_directory(tmp_path, monkeypatch):
    folder = _inputs(tmp_path)
    monkeypatch.chdir(folder)

    result = _run(tmp_path)

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "round1.csv").exists()


# --- annotate -----------------------------------------------------------------------


def _annotate(*args):
    return CliRunner().invoke(cli, ["annotate", "--backend", "ersilia", *args])


def test_annotate_prints_values_priorities_and_the_cutoff_summary(tmp_path):
    folder = _inputs(tmp_path, [QED, ["rdkit", "mw", "0-1000", "range", "secondary"]])

    result = _annotate("--path-to-models", str(folder), "--smiles", "CCO")

    assert result.exit_code == 0, result.output
    assert "qed [primary] (>= 0)" in result.output
    assert "mw [secondary] (0-1000)" in result.output
    assert "1 of 1 primary cutoffs satisfied" in result.output
    assert "guard-rails all satisfied" in result.output


def test_annotate_reports_a_broken_guard_rail(tmp_path):
    folder = _inputs(tmp_path, [QED, ["rdkit", "mw", "250-450", "range", "secondary"]])

    result = _annotate("--path-to-models", str(folder), "--smiles", "CCO")

    assert "guard-rails NOT all satisfied" in result.output


def test_annotate_needs_only_the_annotators(tmp_path):
    folder = tmp_path / "only"
    folder.mkdir()
    _write_csv(folder / "annotators.csv", ANNOTATOR_HEADER, [QED])

    result = _annotate("--path-to-models", str(folder), "--smiles", "CCO")

    assert result.exit_code == 0, result.output


def test_annotate_a_single_in_house_model_prints_the_raw_value():
    result = _annotate("--source", "rdkit", "--model", "mw", "--smiles", "CCO")

    assert result.exit_code == 0, result.output
    assert "mw: 46.069" in result.output
    assert "cutoff" not in result.output


def test_annotate_takes_the_annotators_as_a_csv_flag(tmp_path):
    folder = _inputs(tmp_path)

    result = _annotate(
        "--annotators", str(folder / "annotators.csv"), "--smiles", "CCO"
    )

    assert result.exit_code == 0, result.output
    assert "qed [primary]" in result.output


def test_annotate_finds_annotators_in_the_working_directory(tmp_path, monkeypatch):
    folder = _inputs(tmp_path)
    monkeypatch.chdir(folder)

    result = _annotate("--smiles", "CCO")

    assert result.exit_code == 0, result.output
    assert "qed [primary]" in result.output


def test_annotate_requires_smiles(tmp_path):
    folder = _inputs(tmp_path)

    result = _annotate("--path-to-models", str(folder))

    assert result.exit_code != 0
    assert "--smiles" in result.output


def test_annotate_reports_errors_in_red_and_exits_nonzero(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()

    result = _annotate("--path-to-models", str(empty), "--smiles", "CCO")

    assert result.exit_code == 1
    assert "annotators is required" in result.output


@pytest.mark.parametrize(
    "args, message",
    [
        (["--source", "rdkit"], "go together"),
        (["--model", "mw"], "go together"),
        (["--source", "rdkit", "--model", "mw", "--config", "c.json"], "not both"),
        (["--source", "rdkit", "--model", "mw", "--annotators", "c.json"], "not both"),
        ([], "annotators is required"),
    ],
)
def test_annotate_rejects_inconsistent_options(tmp_path, monkeypatch, args, message):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "c.json").write_text("{}")

    result = _annotate("--smiles", "CCO", *args)

    assert result.exit_code == 1
    assert message in result.output


# --- update-catalog --------------------------------------------------------------------------

GENERATOR = catalog.CatalogEntry(
    "generators", "eos9taz", "MoLeR", "Sampling", "Generation"
)
ANNOTATOR = catalog.CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", "")


def _catalog(*entries, built="2026-10-01"):
    return catalog.Catalog(built, "0.1.62", tuple(entries))


def _update_catalog(monkeypatch, outcome):
    def update(*args, **kwargs):
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(catalog, "update_catalog", update)
    return CliRunner().invoke(cli, ["update-catalog"])


def test_update_catalog_reports_an_up_to_date_catalog(monkeypatch):
    current = _catalog(GENERATOR, ANNOTATOR)
    outcome = catalog.CatalogUpdate(
        current, current, catalog.CatalogDiff((), ()), False
    )

    result = _update_catalog(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert (
        "Catalog is up to date: 1 generators, 0 similarity, 1 annotators"
        in result.output
    )
    assert "2026-10-01" in result.output
    assert "Commit" not in result.output


def test_update_catalog_reports_new_and_removed_models(monkeypatch):
    diff = catalog.CatalogDiff(added=(GENERATOR,), removed=(ANNOTATOR,))
    outcome = catalog.CatalogUpdate(
        _catalog(GENERATOR, built="2026-10-29"), _catalog(ANNOTATOR), diff, True
    )

    result = _update_catalog(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert "Catalog updated (was built 2026-10-01)" in result.output
    assert "new generators: eos9taz (MoLeR)" in result.output
    assert "removed annotators: eos4zfy (MAIP)" in result.output
    assert "Commit" in result.output


def test_update_catalog_creates_a_missing_catalog(monkeypatch):
    diff = catalog.CatalogDiff(added=(GENERATOR,), removed=())
    outcome = catalog.CatalogUpdate(_catalog(GENERATOR), None, diff, True)

    result = _update_catalog(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert "Catalog created" in result.output


def test_update_catalog_fails_cleanly_when_the_hub_cannot_be_reached(monkeypatch):
    result = _update_catalog(
        monkeypatch, RuntimeError("could not query the Ersilia Model Hub")
    )

    assert result.exit_code == 1
    assert "existing catalog was left untouched" in result.output


def test_update_catalog_fails_cleanly_when_the_catalog_cannot_be_written(monkeypatch):
    result = _update_catalog(monkeypatch, PermissionError("read-only file system"))

    assert result.exit_code == 1
    assert (
        "could not write" in result.output and "development checkout" in result.output
    )


def test_update_catalog_writes_the_catalog_file_end_to_end(tmp_path, monkeypatch):
    path = tmp_path / "catalog.csv"
    monkeypatch.setattr(catalog, "CATALOG_PATH", path)

    def hub(*args):
        if args[0] == "--more":
            kind = {"sampling": "Sampling", "annotation": "Annotation"}[args[2]]
            ident = {"sampling": "eos9taz", "annotation": "eos4zfy"}[args[2]]
            return [{"Identifier": ident, "Title": "T", "Task": kind}]
        return {"Subtask": "Generation"}

    monkeypatch.setattr(catalog, "_run_ersilia", hub)

    result = CliRunner().invoke(cli, ["update-catalog"])
    again = CliRunner().invoke(cli, ["update-catalog"])

    assert result.exit_code == 0 and "Catalog created: 1 generators" in result.output
    assert [e.identifier for e in catalog.read_catalog(path).entries] == [
        "eos9taz",
        "eos4zfy",
    ]
    assert "up to date" in again.output


# --- a run inside a workspace ---------------------------------------------------------


def _workspace(tmp_path, built=("eos9taz", "eos8nyg")):
    """A workspace whose models are cloned and whose environments are built."""
    create_workspace(tmp_path)
    _write_csv(tmp_path / "input" / "annotators.csv", ANNOTATOR_HEADER, [QED])
    for model in built:
        (tmp_path / "models" / model).mkdir(parents=True)
        env = tmp_path / "envs" / model
        env.mkdir(parents=True)
        write_marker(env, model=model)
    return tmp_path


def test_a_run_in_a_workspace_writes_a_dated_round_folder(tmp_path, monkeypatch):
    root = _workspace(tmp_path)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(
        cli, ["run", "--mode", "sequential", "--seed-smiles", "CCO", "--n-rounds", "1"]
    )

    assert result.exit_code == 0, result.output
    rounds = list((root / "rounds").iterdir())
    assert len(rounds) == 1
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", rounds[0].name)
    assert (rounds[0] / "summary.csv").exists() and (rounds[0] / "round1.csv").exists()


def test_a_run_saves_what_it_resolved_next_to_its_results(tmp_path, monkeypatch):
    root = _workspace(tmp_path)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(
        cli, ["run", "--mode", "sequential", "--seed-smiles", "CCO", "--n-rounds", "1"]
    )

    assert result.exit_code == 0, result.output
    [round_dir] = list((root / "rounds").iterdir())
    saved = json.loads((round_dir / "inputs.json").read_text())
    assert saved["similarity"] == "eos8nyg" and saved["generators"] == ["eos9taz"]
    assert saved["annotators"] == [
        {
            "model": "qed",
            "cutoff": 0.0,
            "direction": "higher",
            "priority": "primary",
            "weight": 1.0,
        }
    ]
    assert saved["path_to_envs"] == str(root / "envs")


def test_two_runs_in_a_row_get_separate_folders(tmp_path, monkeypatch):
    root = _workspace(tmp_path)
    monkeypatch.chdir(root)
    command = ["run", "--mode", "sequential", "--seed-smiles", "CCO", "--n-rounds", "1"]

    first = CliRunner().invoke(cli, command)
    second = CliRunner().invoke(cli, command)

    assert first.exit_code == 0 and second.exit_code == 0, second.output
    assert len(list((root / "rounds").iterdir())) == 2


def test_an_explicit_output_dir_still_wins_in_a_workspace(tmp_path, monkeypatch):
    root = _workspace(tmp_path)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--mode",
            "sequential",
            "--seed-smiles",
            "CCO",
            "--n-rounds",
            "1",
            "--output-dir",
            str(tmp_path / "elsewhere"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "elsewhere" / "summary.csv").exists()
    assert not list((root / "rounds").iterdir())


def test_outside_a_workspace_the_output_dir_is_still_required(tmp_path, monkeypatch):
    folder = _inputs(tmp_path)
    monkeypatch.chdir(tmp_path)

    result = _run_without_output(
        "--path-to-models", str(folder), "--backend", "ersilia"
    )

    assert result.exit_code == 1
    assert "--output-dir is required outside a workspace" in result.output


def test_a_model_without_a_built_environment_says_to_run_fetch(tmp_path, monkeypatch):
    root = _workspace(tmp_path, built=["eos9taz"])  # eos8nyg has neither clone nor env
    monkeypatch.chdir(root)

    result = CliRunner().invoke(
        cli, ["run", "--mode", "sequential", "--seed-smiles", "CCO"]
    )

    assert result.exit_code == 1
    assert "chemsampler fetch" in result.output
