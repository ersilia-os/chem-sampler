import csv
import importlib.metadata
import json
import re

import pandas as pd
import pytest
from click.testing import CliRunner

from chemsampler import config
from chemsampler.cli import cli
from chemsampler.hub import catalog
from chemsampler.hub.envs import write_marker
from chemsampler.project import create_project

from .conftest import make_rdkit_env


class FakeHubModel:
    """Stands in for HubModel so fake Hub models need no real Ersilia calls."""

    def __init__(self, model_id, models_dir=None, envs_dir=None):
        self.model_id = model_id
        self.elapsed = 0.0

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


ANNOTATOR_HEADER = ["source", "model", "version", "cutoff", "direction", "priority"]
QED = ["rdkit", "qed", "", "0.0", "higher", "primary"]


def _fetch(models_dir, envs_dir, *model_ids):
    """Make models look fetched: a clone folder and a finished environment each."""
    for model in model_ids:
        (models_dir / model).mkdir(parents=True)
        env = envs_dir / model
        env.mkdir(parents=True)
        write_marker(env, model=model)


def _inputs(tmp_path, annotators=(QED,), fetched=True):
    """
    A project with these inputs and, unless `fetched` is False, its two Hub
    models and the rdkit environment fetched into its own folders.
    """
    root = tmp_path / "proj"
    create_project(root, tmp_path / "models", tmp_path / "envs")
    folder = root / "input"
    _write_csv(folder / "similarity.csv", ["model", "version"], [["eos8nyg", "v1.0.0"]])
    _write_csv(folder / "generators.csv", ["model", "version"], [["eos9taz", "v1.0.0"]])
    _write_csv(folder / "annotators.csv", ANNOTATOR_HEADER, list(annotators))
    if fetched:
        _fetch(tmp_path / "models", tmp_path / "envs", "eos8nyg", "eos9taz")
        make_rdkit_env(tmp_path / "envs")
    return root


def _run_dir(root):
    """The one dated folder a run wrote under the project's runs/."""
    [folder] = list((root / "runs").iterdir())
    return folder


def _seed_dir(root, k=1):
    """The folder of the run's k-th seed, with its rounds and summary."""
    return _run_dir(root) / f"seed{k}"


def _run(root, *args):
    return CliRunner().invoke(
        cli,
        [
            "run",
            "-p",
            str(root),
            "--seed-smiles",
            "CCO",
            "--n-rounds",
            "1",
            *args,
        ],
    )


# --- run ----------------------------------------------------------------------------


def test_run_writes_results_and_merges_generator_and_similarity_candidates(tmp_path):
    root = _inputs(tmp_path)

    result = _run(root)

    assert result.exit_code == 0, result.output
    assert "Best candidate" in result.output
    assert (_seed_dir(root) / "summary.csv").exists()
    round1 = pd.read_csv(_seed_dir(root) / "round1.csv")
    assert round1["source"].tolist() == ["eos9taz|eos8nyg"]


def test_a_guard_rail_annotator_acts_as_one_through_the_cli(tmp_path):
    root = _inputs(tmp_path, [QED, ["rdkit", "mw", "", "0-10", "range", "guard-rail"]])

    result = _run(root)

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(_seed_dir(root) / "round1.csv")
    assert round1["guardrails_satisfied"].tolist() == [False]  # CCN weighs about 45 Da


def test_one_primary_runs_sequentially_without_a_mode(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "mw", "", "0-1000", "range", "guard-rail"]]
    )

    result = _run(root)

    assert result.exit_code == 0, result.output
    saved = json.loads((_run_dir(root) / "inputs.json").read_text())
    assert saved["run"]["mode"] == "sequential"


@pytest.mark.parametrize("mode", ["joint", "sequential"])
def test_one_primary_takes_no_mode(tmp_path, mode):
    root = _inputs(tmp_path)

    result = _run(root, "--mode", mode)

    assert result.exit_code == 1
    assert "a run is always sequential: drop --mode" in result.output
    assert not list((root / "runs").iterdir())


def test_several_primaries_need_a_mode(tmp_path):
    root = _inputs(tmp_path, [QED, ["rdkit", "mw", "", "0", "higher", "primary"]])

    refused = _run(root)
    weighted = _run(root, "--mode", "weighted")

    assert refused.exit_code == 1
    assert "with 2 primary annotators a mode is needed (--mode)" in refused.output
    assert weighted.exit_code == 0, weighted.output
    saved = json.loads((_run_dir(root) / "inputs.json").read_text())
    assert saved["run"]["mode"] == "weighted"


def test_run_needs_the_inputs_and_reports_errors_in_red(tmp_path):
    root = _inputs(tmp_path)
    (root / "input" / "similarity.csv").unlink()
    (root / "input" / "generators.csv").unlink()

    result = _run(root)

    assert result.exit_code == 1
    assert "give at least one of similarity and generators" in result.output
    assert not list((root / "runs").iterdir())


def test_a_run_with_nothing_fetched_says_to_run_fetch(tmp_path):
    root = _inputs(tmp_path, fetched=False)

    result = CliRunner().invoke(cli, ["run", "-p", str(root), "--seed-smiles", "CCO"])

    assert result.exit_code == 1
    assert "no models fetched yet" in result.output
    assert "chemsampler fetch" in result.output


def test_a_run_needs_a_seed(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "tanimoto", "", "0.4", "higher", "guard-rail"]]
    )

    result = CliRunner().invoke(cli, ["run", "-p", str(root)])

    assert result.exit_code == 1
    assert "give a seed: input/smiles.csv" in result.output
    assert "--seed-smiles" in result.output
    assert not list((root / "runs").iterdir())


@pytest.mark.parametrize("smiles", ["CCO,CCN", "CCO CCN", "not-a-molecule"])
def test_seed_smiles_takes_one_smiles(tmp_path, smiles):
    root = _inputs(tmp_path)

    result = CliRunner().invoke(cli, ["run", "-p", str(root), "--seed-smiles", smiles])

    assert result.exit_code == 1
    assert "--seed-smiles takes one SMILES" in result.output
    assert "put several seeds in input/smiles.csv" in result.output
    assert not list((root / "runs").iterdir())


def test_a_run_records_how_it_was_launched(tmp_path):
    root = _inputs(tmp_path)

    result = _run(root, "--tolerance", "0.1")

    assert result.exit_code == 0, result.output
    saved = json.loads((_run_dir(root) / "inputs.json").read_text())
    # the models are still there, with their versions
    assert saved["similarity"] == [{"model": "eos8nyg", "version": "v1.0.0"}]
    assert saved["run"] == {
        "chemsampler_version": importlib.metadata.version("chemsampler"),
        "mode": "sequential",
        "seeds": ["CCO"],
        "n_rounds": 1,
        "n_steps": None,
        "tolerance": 0.1,
        "hpc": None,
        "verbose": False,
        "status": "finished",
    }


def test_an_interrupted_run_keeps_its_record_with_the_partial_results(
    tmp_path, monkeypatch
):
    from chemsampler import optimize

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(optimize, "hill_climb", interrupted)
    root = _inputs(tmp_path)

    result = _run(root)

    assert result.exit_code == 0, result.output
    assert "Partial results written" in result.output
    saved = json.loads((_run_dir(root) / "inputs.json").read_text())
    assert saved["run"]["mode"] == "sequential"
    assert (_seed_dir(root) / "summary.csv").exists()


#: CCN (45.1 Da), the fake models' only candidate, beats the seed CCO (46.1 Da),
#: and the cutoff is out of reach, so the run goes on to round 2.
MW_LOWER = ["rdkit", "mw", "", "10", "lower", "primary"]


def _stop_in_round_2(monkeypatch, error):
    """Make the fake models raise `error` once round 1 has run both of them."""
    calls = []

    def run(self, smiles_list):
        calls.append(self.model_id)
        if len(calls) > 2:
            raise error
        return pd.DataFrame({"smi_0": ["CCN"]})

    monkeypatch.setattr(FakeHubModel, "run", run)


def test_ctrl_c_keeps_the_rounds_finished_before_it(tmp_path, monkeypatch):
    _stop_in_round_2(monkeypatch, KeyboardInterrupt())
    root = _inputs(tmp_path, [MW_LOWER])

    result = _run(root, "--n-rounds", "3")

    assert result.exit_code == 0, result.output
    assert "Partial results written" in result.output
    summary = pd.read_csv(_seed_dir(root) / "summary.csv")
    assert summary["round"].tolist() == [0, 1]
    assert (_seed_dir(root) / "round1.csv").exists()


def test_a_failing_model_keeps_the_rounds_finished_before_it(tmp_path, monkeypatch):
    _stop_in_round_2(monkeypatch, RuntimeError("eos9taz: run.sh failed (exit 1)"))
    root = _inputs(tmp_path, [MW_LOWER])

    result = _run(root, "--n-rounds", "3")

    assert result.exit_code == 1
    assert "eos9taz: run.sh failed (exit 1)" in result.output
    assert f"Finished rounds written to {_run_dir(root)}" in result.output
    summary = pd.read_csv(_seed_dir(root) / "summary.csv")
    assert summary["round"].tolist() == [0, 1]


# --- resume -------------------------------------------------------------------------

#: Heavier molecules each round, so every round improves on the last.
MW_HIGHER = ["rdkit", "mw", "", "1000", "higher", "primary"]
HEAVIER = ["CCCC", "CCCCC", "CCCCCC"]


def _script(monkeypatch, *steps):
    """The generator's answers, one per call: a SMILES, an error, or a function."""
    steps = list(steps)

    def run(self, smiles_list):
        step = steps.pop(0)
        if isinstance(step, BaseException):
            raise step
        if callable(step):
            step()
        return pd.DataFrame({"smi_0": [step]})

    monkeypatch.setattr(FakeHubModel, "run", run)
    return steps


def _generator_project(tmp_path, annotators=(MW_HIGHER,)):
    root = _inputs(tmp_path, list(annotators))
    (root / "input" / "similarity.csv").unlink()  # one model call per round
    return root


def _three_rounds(root):
    return _run(root, "--n-rounds", "3")


def _resume(root, *args):
    return CliRunner().invoke(cli, ["resume", "-p", str(root), *args])


def _status(root):
    return json.loads((_run_dir(root) / "inputs.json").read_text())["run"]["status"]


def _same_results(here, there):
    """Two runs' folders hold the same tables, seed by seed, to the last digit."""
    here, there = here / "seed1", there / "seed1"
    names = sorted(p.name for p in here.glob("*.csv"))
    assert names == sorted(p.name for p in there.glob("*.csv"))
    for name in names:
        pd.testing.assert_frame_equal(
            pd.read_csv(here / name, float_precision="round_trip"),
            pd.read_csv(there / name, float_precision="round_trip"),
            check_exact=True,
        )


def test_a_stopped_run_resumes_to_what_an_uninterrupted_run_gives(
    tmp_path, monkeypatch
):
    whole = _generator_project(tmp_path / "whole")
    _script(monkeypatch, *HEAVIER)
    assert _three_rounds(whole).exit_code == 0
    root = _generator_project(tmp_path / "stopped")
    _script(monkeypatch, *HEAVIER[:2], KeyboardInterrupt())

    stopped = _three_rounds(root)

    assert stopped.exit_code == 0, stopped.output
    name = _run_dir(root).name
    assert f"Continue it with `chemsampler resume {name}`" in stopped.output
    assert _status(root) == "stopped"
    assert not (_run_dir(root) / "run.lock").exists()
    left = _script(monkeypatch, HEAVIER[2])

    resumed = _resume(root)

    assert resumed.exit_code == 0, resumed.output
    assert f"Resuming {name} (stopped)." in resumed.output
    assert left == []  # only round 3 ran a model
    assert _status(root) == "finished"
    _same_results(_run_dir(whole), _run_dir(root))


def test_resume_uses_the_inputs_the_run_started_with(tmp_path, monkeypatch):
    root = _generator_project(tmp_path)
    _script(monkeypatch, HEAVIER[0], KeyboardInterrupt())
    _three_rounds(root)
    _write_csv(root / "input" / "annotators.csv", ANNOTATOR_HEADER, [QED])
    _script(monkeypatch, *HEAVIER[1:])

    result = _resume(root, _run_dir(root).name)

    assert result.exit_code == 0, result.output
    assert "mw" in pd.read_csv(_seed_dir(root) / "round3.csv").columns


def test_sigterm_stops_a_run_as_ctrl_c_does(tmp_path, monkeypatch):
    import os
    import signal

    root = _generator_project(tmp_path)
    _script(monkeypatch, HEAVIER[0], lambda: os.kill(os.getpid(), signal.SIGTERM))

    result = _three_rounds(root)

    assert result.exit_code == 0, result.output
    assert "Partial results written" in result.output
    assert _status(root) == "stopped"
    assert pd.read_csv(_seed_dir(root) / "summary.csv")["round"].tolist() == [0, 1]


def test_a_failed_run_can_be_resumed_once_fixed(tmp_path, monkeypatch):
    root = _generator_project(tmp_path)
    _script(monkeypatch, HEAVIER[0], RuntimeError("eos9taz: out of memory"))

    failed = _three_rounds(root)

    assert failed.exit_code == 1
    assert "Once fixed, continue it with `chemsampler resume" in failed.output
    assert _status(root) == "failed"
    _script(monkeypatch, *HEAVIER[1:])
    assert _resume(root).exit_code == 0
    assert _status(root) == "finished"


def test_a_finished_run_is_not_run_again(tmp_path, monkeypatch):
    root = _generator_project(tmp_path)
    _script(monkeypatch, *HEAVIER)
    _three_rounds(root)
    _script(monkeypatch)  # any model call would fail

    result = _resume(root)

    assert result.exit_code == 0, result.output
    assert "already finished" in result.output


def test_resume_explains_what_it_cannot_resume(tmp_path, monkeypatch):
    root = _generator_project(tmp_path)
    assert "no runs in" in _resume(root).output
    old = root / "runs" / "2026-01-01_00-00-00"
    old.mkdir()
    (old / "inputs.json").write_text(json.dumps({"run": {"mode": "weighted"}}))

    result = _resume(root)

    assert result.exit_code == 1
    assert "made before chemsampler could resume runs" in result.output


def test_resume_refuses_a_run_still_going_here(tmp_path, monkeypatch):
    import os
    import socket

    root = _generator_project(tmp_path)
    _script(monkeypatch, HEAVIER[0], KeyboardInterrupt())
    _three_rounds(root)
    lock = _run_dir(root) / "run.lock"
    held = {"host": socket.gethostname(), "pid": os.getpid(), "job": None}
    lock.write_text(json.dumps(held))

    result = _resume(root)

    assert result.exit_code == 1
    assert "still going on this machine" in result.output
    assert json.loads(lock.read_text()) == held


def test_run_can_use_a_similarity_model_without_generators(tmp_path):
    root = _inputs(tmp_path)
    (root / "input" / "generators.csv").unlink()

    result = _run(root)

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(_seed_dir(root) / "round1.csv")
    assert round1["source"].tolist() == ["eos8nyg"]


def test_run_can_use_generators_without_a_similarity_model(tmp_path):
    root = _inputs(tmp_path)
    (root / "input" / "similarity.csv").unlink()

    result = _run(root)

    assert result.exit_code == 0, result.output
    round1 = pd.read_csv(_seed_dir(root) / "round1.csv")
    assert round1["source"].tolist() == ["eos9taz"]


# --- annotate -----------------------------------------------------------------------


# --- the seed (input/smiles.csv) and the Tanimoto annotator ---------------------------

TANIMOTO = ["rdkit", "tanimoto", "", "0.5", "higher", "guard-rail"]


def _seed(root, text):
    (root / "input" / "smiles.csv").write_text(text)


def _run_from_csv(root, *args):
    """A run that takes its seeds from input/smiles.csv, if anything."""
    return CliRunner().invoke(cli, ["run", "-p", str(root), "--n-rounds", "1", *args])


def _record(root):
    return json.loads((_run_dir(root) / "inputs.json").read_text())


def test_the_seed_comes_from_smiles_csv(tmp_path):
    root = _inputs(tmp_path)
    _seed(root, "smiles\nCCO\n")

    result = _run_from_csv(root)

    assert result.exit_code == 0, result.output
    assert _record(root)["run"]["seeds"] == ["CCO"]
    assert (_run_dir(root) / "input" / "smiles.csv").read_text() == "smiles\nCCO\n"
    assert pd.read_csv(_seed_dir(root) / "round0.csv")["smiles"].tolist() == ["CCO"]


def test_seed_smiles_wins_over_smiles_csv(tmp_path):
    root = _inputs(tmp_path)
    _seed(root, "smiles\nCCC\n")

    result = _run(root)  # --seed-smiles CCO

    assert result.exit_code == 0, result.output
    assert "--seed-smiles given: ignoring input/smiles.csv" in result.output
    assert _record(root)["run"]["seeds"] == ["CCO"]


@pytest.mark.parametrize(
    "text, message",
    [
        ("smiles\nCCO\nCCC\nOCC\n", "rows 1 and 3 are the same molecule"),
        ("smiles\n", "smiles.csv has no molecule"),
        ("smiles\nCCO\nnot-a-molecule\n", "row 2: RDKit cannot read the SMILES"),
        ("smiles\nCCO CCN\n", "row 1: RDKit cannot read the SMILES 'CCO CCN'"),
        ("seed\nCCO\n", "smiles.csv"),
    ],
)
def test_a_bad_smiles_csv_is_refused(tmp_path, text, message):
    root = _inputs(tmp_path)
    _seed(root, text)

    result = _run_from_csv(root)

    assert result.exit_code == 1
    assert message in result.output
    assert not (root / "runs").exists() or not any((root / "runs").iterdir())


def test_a_tanimoto_guard_rail_through_the_cli(tmp_path):
    root = _inputs(tmp_path, [QED, TANIMOTO])
    _seed(root, "smiles\nCCO\n")

    result = _run_from_csv(root)

    assert result.exit_code == 0, result.output
    round0 = pd.read_csv(_seed_dir(root) / "round0.csv")
    assert round0["tanimoto"].tolist() == [1.0]
    round1 = pd.read_csv(_seed_dir(root) / "round1.csv")  # CCN, 0.333 to CCO
    assert round1["tanimoto"].tolist() == [round1["tanimoto_to_seed"][0]]
    assert round1["guardrails_satisfied"].tolist() == [False]


def test_each_seed_is_climbed_in_its_own_folder_with_tanimoto_to_its_own_seed(
    tmp_path,
):
    root = _inputs(tmp_path, [QED, TANIMOTO])
    _seed(root, "smiles\nCCO\nc1ccccc1\n")

    result = _run_from_csv(root)

    assert result.exit_code == 0, result.output
    assert "Seed 1 of 2: CCO (seed1/)" in result.output
    assert "Seed 2 of 2: c1ccccc1 (seed2/)" in result.output
    assert "Best candidate from CCO (seed1/)" in result.output
    assert "Best candidate from c1ccccc1 (seed2/)" in result.output
    assert _record(root)["run"]["seeds"] == ["CCO", "c1ccccc1"]
    assert sorted(p.name for p in _run_dir(root).glob("seed*")) == ["seed1", "seed2"]
    for k, seed in [(1, "CCO"), (2, "c1ccccc1")]:
        round0 = pd.read_csv(_seed_dir(root, k) / "round0.csv")
        assert round0["smiles"].tolist() == [seed] and round0["tanimoto"][0] == 1.0
        assert (_seed_dir(root, k) / "summary.csv").is_file()
    # the same candidate (CCN), compared with each seed
    near, far = (pd.read_csv(_seed_dir(root, k) / "round1.csv") for k in (1, 2))
    assert near["tanimoto"][0] > far["tanimoto"][0]
    assert near["tanimoto"][0] == near["tanimoto_to_seed"][0]


def test_a_run_stopped_in_its_second_seed_resumes_there(tmp_path, monkeypatch):
    root = _generator_project(tmp_path)
    _seed(root, "smiles\nCCO\nCCN\n")
    # seed 1: rounds 1-3; seed 2: round 1, then Ctrl-C in round 2
    _script(monkeypatch, *HEAVIER, "CCCCCCC", KeyboardInterrupt())

    stopped = _run_from_csv(root, "--n-rounds", "3")

    assert stopped.exit_code == 0, stopped.output
    assert f"Partial results written to {_seed_dir(root, 2)}" in stopped.output
    assert _status(root) == "stopped"
    seed1 = pd.read_csv(_seed_dir(root, 1) / "summary.csv")
    assert seed1["round"].tolist() == [0, 1, 2, 3]
    assert pd.read_csv(_seed_dir(root, 2) / "summary.csv")["round"].tolist() == [0, 1]
    before = (_seed_dir(root, 1) / "round3.csv").read_text()
    calls = _script(monkeypatch, "CCCCCCCC", "CCCCCCCCC")

    resumed = _resume(root)

    assert resumed.exit_code == 0, resumed.output
    assert calls == []  # only seed 2's rounds 2 and 3 called the generator
    assert _status(root) == "finished"
    assert (_seed_dir(root, 1) / "round3.csv").read_text() == before
    seed2 = pd.read_csv(_seed_dir(root, 2) / "summary.csv")
    assert seed2["round"].tolist() == [0, 1, 2, 3]


def test_resume_refuses_a_run_made_before_seed_folders(tmp_path):
    root = _inputs(tmp_path)
    assert _run(root).exit_code == 0
    record = _record(root)
    record["run"]["status"] = "stopped"
    record["run"]["seed_smiles"] = record["run"].pop("seeds")[0]
    (_run_dir(root) / "inputs.json").write_text(json.dumps(record))

    result = _resume(root)

    assert result.exit_code == 1
    assert "before each seed had its own folder: start it again" in result.output


@pytest.mark.parametrize(
    "option", ["--tanimoto-cutoff", "--tanimoto-direction", "--original-seed-smiles"]
)
def test_the_tanimoto_options_are_gone(tmp_path, option):
    result = _run(_inputs(tmp_path), option, "0.5")

    assert result.exit_code == 2
    assert "No such option" in result.output


def test_a_stopped_run_with_a_tanimoto_annotator_resumes(tmp_path, monkeypatch):
    tanimoto = ["rdkit", "tanimoto", "", "0", "higher", "guard-rail"]
    whole = _generator_project(tmp_path / "whole", [MW_HIGHER, tanimoto])
    _script(monkeypatch, *HEAVIER)
    assert _three_rounds(whole).exit_code == 0
    root = _generator_project(tmp_path / "stopped", [MW_HIGHER, tanimoto])
    _script(monkeypatch, *HEAVIER[:2], KeyboardInterrupt())
    assert _three_rounds(root).exit_code == 0
    _script(monkeypatch, HEAVIER[2])

    resumed = _resume(root)

    assert resumed.exit_code == 0, resumed.output
    _same_results(_run_dir(whole), _run_dir(root))


@pytest.mark.parametrize(
    "key, value", [("tanimoto_cutoff", 0.3), ("original_seed_smiles", "CCO")]
)
def test_resume_refuses_a_run_made_with_the_removed_options(
    tmp_path, monkeypatch, key, value
):
    root = _generator_project(tmp_path)
    _script(monkeypatch, HEAVIER[0], KeyboardInterrupt())
    _three_rounds(root)
    record = _record(root)
    record["run"][key] = value
    (_run_dir(root) / "inputs.json").write_text(json.dumps(record))

    result = _resume(root)

    assert result.exit_code == 1
    assert "made by an older chemsampler" in result.output


def test_annotate_measures_tanimoto_against_the_seed(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "tanimoto", "", "0.2", "higher", "guard-rail"]]
    )
    _seed(root, "smiles\nCCO\n")

    result = _annotate("-p", str(root), "--smiles", "CO")

    assert result.exit_code == 0, result.output
    assert "tanimoto [guard-rail] (>= 0.2): 0.2857142857142857" in result.output
    assert "guard-rails all satisfied" in result.output


def test_annotate_without_a_seed_leaves_tanimoto_out(tmp_path):
    root = _inputs(tmp_path, [QED, TANIMOTO])

    result = _annotate("-p", str(root), "--smiles", "CO")

    assert result.exit_code == 0, result.output
    assert "tanimoto [guard-rail] (>= 0.5): needs a seed" in result.output
    assert "1 of 1 primary cutoffs satisfied" in result.output


def test_annotate_a_single_tanimoto_uses_the_project_seed(tmp_path):
    root = _inputs(tmp_path)
    args = ["-p", str(root), "--source", "rdkit", "--model", "tanimoto"]

    without = _annotate(*args, "--smiles", "CO")
    _seed(root, "smiles\nCCO\n")
    seeded = _annotate(*args, "--smiles", "CO")

    assert without.exit_code == 1
    assert "tanimoto compares candidates with the seed" in without.output
    assert seeded.exit_code == 0, seeded.output
    assert "tanimoto: 0.2857142857142857" in seeded.output


def test_annotate_with_several_seeds_shows_tanimoto_per_seed_uncounted(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "tanimoto", "", "0.2", "higher", "guard-rail"]]
    )
    _seed(root, "smiles\nCCO\nCO\n")
    single = ["-p", str(root), "--source", "rdkit", "--model", "tanimoto"]

    result = _annotate("-p", str(root), "--smiles", "CO")
    one_model = _annotate(*single, "--smiles", "CO")

    assert result.exit_code == 0, result.output
    assert (
        "tanimoto [guard-rail] (>= 0.2): tanimoto vs seed 1: 0.2857142857142857; "
        "tanimoto vs seed 2: 1.0; not counted"
    ) in result.output
    assert "qed [primary]" in result.output
    assert "1 of 1 primary cutoffs satisfied" in result.output
    assert "guard-rails" not in result.output  # none counted
    assert one_model.exit_code == 0, one_model.output
    assert "tanimoto vs seed 1: 0.2857142857142857" in one_model.output
    assert "tanimoto vs seed 2: 1.0" in one_model.output


def _annotate(*args):
    return CliRunner().invoke(cli, ["annotate", *args])


def test_annotate_prints_values_priorities_and_the_cutoff_summary(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "mw", "", "0-1000", "range", "guard-rail"]]
    )

    result = _annotate("-p", str(root), "--smiles", "CCO")

    assert result.exit_code == 0, result.output
    assert "qed [primary] (>= 0)" in result.output
    assert "mw [guard-rail] (0-1000)" in result.output
    assert "1 of 1 primary cutoffs satisfied" in result.output
    assert "guard-rails all satisfied" in result.output


def test_annotate_reports_a_broken_guard_rail(tmp_path):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "mw", "", "250-450", "range", "guard-rail"]]
    )

    result = _annotate("-p", str(root), "--smiles", "CCO")

    assert "guard-rails NOT all satisfied" in result.output


def test_annotate_needs_only_the_annotators(tmp_path):
    root = _inputs(tmp_path)
    (root / "input" / "similarity.csv").unlink()
    (root / "input" / "generators.csv").unlink()

    result = _annotate("-p", str(root), "--smiles", "CCO")

    assert result.exit_code == 0, result.output


def test_annotate_a_single_in_house_model_needs_no_project(
    tmp_path, store_root, monkeypatch
):
    monkeypatch.chdir(tmp_path)

    missing = _annotate("--source", "rdkit", "--model", "mw", "--smiles", "CCO")
    make_rdkit_env(store_root / "envs")  # the store's, as no project names one
    result = _annotate("--source", "rdkit", "--model", "mw", "--smiles", "CCO")

    assert missing.exit_code == 1 and "`chemsampler fetch rdkit`" in missing.output
    assert result.exit_code == 0, result.output
    assert "mw: 46.069" in result.output
    assert "cutoff" not in result.output


def test_annotate_a_single_model_looks_in_the_project_given(tmp_path, monkeypatch):
    root = tmp_path / "proj"
    create_project(root, models=tmp_path / "mine", envs=tmp_path / "mine-envs")
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(
        cli,
        ["annotate", "-p", str(root), "--source", "ersilia", "--model", "eos4zfy"]
        + ["--smiles", "CCO"],
    )

    # the project's folder, not the store: nothing is fetched there yet
    assert result.exit_code == 1
    assert f"{tmp_path / 'mine'} does not exist" in result.output


def test_annotate_a_single_model_that_was_not_fetched_says_to_fetch_it(tmp_path):
    root = _inputs(tmp_path)  # eos8nyg and eos9taz are fetched, eos4zfy is not

    result = CliRunner().invoke(
        cli,
        ["annotate", "-p", str(root), "--source", "ersilia", "--model", "eos4zfy"]
        + ["--smiles", "CCO"],
    )

    assert result.exit_code == 1
    assert "no folder for ['eos4zfy']" in result.output
    assert "chemsampler fetch" in result.output


def test_annotate_requires_smiles(tmp_path):
    root = _inputs(tmp_path)

    result = _annotate("-p", str(root))

    assert result.exit_code != 0
    assert "--smiles" in result.output


def test_annotate_reports_errors_in_red_and_exits_nonzero(tmp_path):
    root = _inputs(tmp_path)
    (root / "input" / "annotators.csv").unlink()

    result = _annotate("-p", str(root), "--smiles", "CCO")

    assert result.exit_code == 1
    assert f"annotators is required: write {root / 'input' / 'annotators.csv'}" in (
        result.output
    )


@pytest.mark.parametrize(
    "args, message",
    [
        (["--source", "rdkit"], "go together"),
        (["--model", "mw"], "go together"),
        ([], "not inside a project"),
    ],
)
def test_annotate_rejects_inconsistent_options(tmp_path, monkeypatch, args, message):
    monkeypatch.chdir(tmp_path)

    result = _annotate("--smiles", "CCO", *args)

    assert result.exit_code == 1
    assert message in result.output


# --- catalog --update --------------------------------------------------------------------------

GENERATOR = catalog.CatalogEntry(
    "generators", "eos9taz", "MoLeR", "Sampling", "Generation"
)
ANNOTATOR = catalog.CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", "")


def _catalog(*entries, built="2026-10-01"):
    return catalog.Catalog(built, "0.1.62", tuple(entries))


def _catalog_update(monkeypatch, outcome):
    def update(*args, **kwargs):
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(catalog, "update_catalog", update)
    return CliRunner().invoke(cli, ["catalog", "--update"])


def test_catalog_update_reports_an_up_to_date_catalog(monkeypatch):
    current = _catalog(GENERATOR, ANNOTATOR)
    outcome = catalog.CatalogUpdate(
        current, current, catalog.CatalogDiff((), ()), False
    )

    result = _catalog_update(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert (
        "Catalog is up to date: 1 generators, 0 similarity, 1 annotators"
        in result.output
    )
    assert "2026-10-01" in result.output
    assert "Commit" not in result.output


def test_catalog_update_reports_new_and_removed_models(monkeypatch):
    diff = catalog.CatalogDiff(added=(GENERATOR,), removed=(ANNOTATOR,))
    outcome = catalog.CatalogUpdate(
        _catalog(GENERATOR, built="2026-10-29"), _catalog(ANNOTATOR), diff, True
    )

    result = _catalog_update(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert "Catalog updated (was built 2026-10-01)" in result.output
    assert "new generators: eos9taz (MoLeR)" in result.output
    assert "removed annotators: eos4zfy (MAIP)" in result.output
    assert "Commit" in result.output


def test_catalog_update_creates_a_missing_catalog(monkeypatch):
    diff = catalog.CatalogDiff(added=(GENERATOR,), removed=())
    outcome = catalog.CatalogUpdate(_catalog(GENERATOR), None, diff, True)

    result = _catalog_update(monkeypatch, outcome)

    assert result.exit_code == 0, result.output
    assert "Catalog created" in result.output


def test_catalog_update_fails_cleanly_when_the_hub_cannot_be_reached(monkeypatch):
    result = _catalog_update(
        monkeypatch, RuntimeError("could not query the Ersilia Model Hub")
    )

    assert result.exit_code == 1
    assert "existing catalog was left untouched" in result.output


def test_catalog_update_fails_cleanly_when_the_catalog_cannot_be_written(monkeypatch):
    result = _catalog_update(monkeypatch, PermissionError("read-only file system"))

    assert result.exit_code == 1
    assert (
        "could not write" in result.output and "development checkout" in result.output
    )


def test_catalog_update_writes_the_catalog_file_end_to_end(tmp_path, monkeypatch):
    path = tmp_path / "catalog.csv"
    monkeypatch.setattr(catalog, "CATALOG_PATH", path)

    def hub(*args):
        if args[0] == "--more":
            kind = {"sampling": "Sampling", "annotation": "Annotation"}[args[2]]
            ident = {"sampling": "eos9taz", "annotation": "eos4zfy"}[args[2]]
            return [{"Identifier": ident, "Title": "T", "Task": kind}]
        return {"Subtask": "Generation"}

    monkeypatch.setattr(catalog, "_run_ersilia", hub)

    result = CliRunner().invoke(cli, ["catalog", "--update"])
    again = CliRunner().invoke(cli, ["catalog", "--update"])

    assert result.exit_code == 0 and "Catalog created: 1 generators" in result.output
    assert [e.identifier for e in catalog.read_catalog(path).entries] == [
        "eos9taz",
        "eos4zfy",
    ]
    assert "up to date" in again.output


# --- a run inside a project ---------------------------------------------------------


def _project(tmp_path, store_root, built=("eos9taz", "eos8nyg")):
    """
    A project, and the machine's store holding the models it uses.

    The project owns the inputs and the results; the clones and the built
    environments are in the store, which is where `fetch` puts them.
    """
    root = tmp_path / "proj"
    create_project(root)
    _write_csv(root / "input" / "annotators.csv", ANNOTATOR_HEADER, [QED])
    _fetch(store_root / "models", store_root / "envs", *built)
    make_rdkit_env(store_root / "envs")
    return root


def test_a_run_in_a_project_writes_a_dated_round_folder(
    tmp_path, store_root, monkeypatch
):
    root = _project(tmp_path, store_root)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(cli, ["run", "--seed-smiles", "CCO", "--n-rounds", "1"])

    assert result.exit_code == 0, result.output
    runs = list((root / "runs").iterdir())
    assert len(runs) == 1
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", runs[0].name)
    assert (runs[0] / "seed1" / "summary.csv").exists() and (
        runs[0] / "seed1" / "round1.csv"
    ).exists()


def test_a_run_saves_what_it_resolved_next_to_its_results(
    tmp_path, store_root, monkeypatch
):
    root = _project(tmp_path, store_root)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(cli, ["run", "--seed-smiles", "CCO", "--n-rounds", "1"])

    assert result.exit_code == 0, result.output
    [round_dir] = list((root / "runs").iterdir())
    saved = json.loads((round_dir / "inputs.json").read_text())
    assert saved["similarity"] == [{"model": "eos8nyg", "version": "live"}]
    assert saved["generators"] == [{"model": "eos9taz", "version": "live"}]
    assert saved["annotators"] == [
        {
            "model": "qed",
            "version": None,
            "cutoff": 0.0,
            "direction": "higher",
            "priority": "primary",
            "weight": 1.0,
        }
    ]
    # which models ran is recorded too: they are the machine's, not the project's
    assert saved["path_to_envs"] == str(store_root / "envs")
    assert saved["path_to_model_repositories"] == str(store_root / "models")


def test_two_runs_in_a_row_get_separate_folders(tmp_path, store_root, monkeypatch):
    root = _project(tmp_path, store_root)
    monkeypatch.chdir(root)
    command = ["run", "--seed-smiles", "CCO", "--n-rounds", "1"]

    first = CliRunner().invoke(cli, command)
    second = CliRunner().invoke(cli, command)

    assert first.exit_code == 0 and second.exit_code == 0, second.output
    assert len(list((root / "runs").iterdir())) == 2


def test_a_model_without_a_built_environment_says_to_run_fetch(
    tmp_path, store_root, monkeypatch
):
    root = _project(
        tmp_path, store_root, built=["eos9taz"]
    )  # eos8nyg: no clone, no env
    monkeypatch.chdir(root)

    result = CliRunner().invoke(cli, ["run", "--seed-smiles", "CCO"])

    assert result.exit_code == 1
    assert "chemsampler fetch" in result.output


# --- fetch ----------------------------------------------------------------------------


@pytest.fixture
def fetched(monkeypatch):
    """Replace what fetch would actually do, recording the plan it was given."""
    from chemsampler.hub import fetch

    seen = {}
    monkeypatch.setattr(fetch, "check_prerequisites", lambda: None)

    def fake_fetch_models(plan, *, update, rebuild, say):
        seen.update(plan=plan, update=update, rebuild=rebuild)
        return [
            fetch.ModelResult(m, clone="cloned", checkpoints="none needed", env="built")
            for m in plan.model_ids
        ]

    monkeypatch.setattr(fetch, "fetch_models", fake_fetch_models)
    return seen


def test_fetch_puts_models_in_the_store_and_says_where(
    tmp_path, store_root, monkeypatch, fetched
):
    root = tmp_path / "proj"
    create_project(root)
    monkeypatch.chdir(root)

    result = CliRunner().invoke(cli, ["fetch", "-f"])

    assert result.exit_code == 0, result.output
    assert fetched["plan"].models_dir == store_root / "models"
    assert fetched["plan"].envs_dir == store_root / "envs"
    # ~/.chemsampler is hidden: the user has to be told where gigabytes went
    assert f"clones:       {store_root / 'models'}" in result.output
    assert f"environments: {store_root / 'envs'}" in result.output
    assert f"logs:         {store_root / 'logs'}" in result.output
    assert "All models are ready." in result.output


def test_fetch_uses_the_folders_the_project_names(tmp_path, monkeypatch, fetched):
    root = tmp_path / "proj"
    create_project(
        root, models=tmp_path / "shared" / "models", envs=tmp_path / "shared" / "envs"
    )
    monkeypatch.chdir(root)

    result = CliRunner().invoke(cli, ["fetch", "-f"])

    assert result.exit_code == 0, result.output
    assert fetched["plan"].models_dir == tmp_path / "shared" / "models"
    assert fetched["plan"].envs_dir == tmp_path / "shared" / "envs"


def test_plain_fetch_asks_before_fetching_everything_listed(
    tmp_path, fetched, monkeypatch
):
    root = _inputs(
        tmp_path, [QED, ["rdkit", "tanimoto", "", "0.4", "higher", "guard-rail"]]
    )

    no = _fetch_ids(root, can_ask=True, answer="n\n", monkeypatch=monkeypatch)
    nobody = _fetch_ids(root, monkeypatch=monkeypatch)
    yes = _fetch_ids(root, can_ask=True, answer="y\n", monkeypatch=monkeypatch)

    listed = (
        "2 model(s): eos9taz v1.0.0, eos8nyg v1.0.0, and the rdkit environment "
        "(qed, tanimoto)"
    )
    assert no.exit_code == 0 and f"input/*.csv list: {listed}" in no.output
    assert "Fetch these? [y/N]: n" in no.output and "Nothing fetched." in no.output
    assert nobody.exit_code == 1 and "run it with -f" in nobody.output
    assert yes.exit_code == 0, yes.output
    assert f"Fetching {listed}" in yes.output
    plan = fetched["plan"]
    assert plan.model_ids == ("eos9taz", "eos8nyg") and plan.rdkit
    assert plan.versions == {"eos9taz": "v1.0.0", "eos8nyg": "v1.0.0"}


def test_plain_fetch_with_f_does_not_ask(tmp_path, fetched, monkeypatch):
    root = _inputs(tmp_path)

    result = _fetch_ids(root, "-f", monkeypatch=monkeypatch)

    assert result.exit_code == 0, result.output
    assert "Fetch these?" not in result.output and fetched["plan"].rdkit


def test_fetch_rdkit_alone_builds_only_its_environment(tmp_path, fetched):
    root = _inputs(tmp_path)
    (root / "input" / "generators.csv").unlink()  # inputs it does not read
    (root / "input" / "similarity.csv").unlink()

    result = _fetch_ids(root, "rdkit")

    assert result.exit_code == 0, result.output
    assert "Fetching the rdkit environment" in result.output
    assert fetched["plan"].model_ids == () and fetched["plan"].rdkit
    assert fetched["plan"].envs_dir == tmp_path / "envs"
    assert sorted(p.name for p in (root / "input").iterdir()) == ["annotators.csv"]


@pytest.mark.parametrize("name", ["mw", "qed", "tanimoto"])
def test_fetch_refuses_rdkit_annotators_as_ids(tmp_path, fetched, name):
    root = _inputs(tmp_path)

    result = _fetch_ids(root, "eos9p57", name, "-f")

    assert result.exit_code == 1
    assert f"{name}: not models but RDKit annotators" in result.output
    assert "use `chemsampler fetch rdkit`" in result.output
    assert "plan" not in fetched


def test_fetch_pulls_and_rebuilds_by_default(tmp_path, fetched):
    root = _inputs(tmp_path)

    result = CliRunner().invoke(cli, ["fetch", "-p", str(root), "-f"])
    old_flag = CliRunner().invoke(cli, ["fetch", "-p", str(root), "--update"])

    assert result.exit_code == 0, result.output
    assert fetched["update"] is True and fetched["rebuild"] is True
    assert old_flag.exit_code == 2 and "No such option" in old_flag.output


# --- fetch with model ids ----------------------------------------------------------


def _fetch_ids(root, *args, can_ask=False, answer=None, monkeypatch=None):
    if monkeypatch is not None:
        from chemsampler.cli import options

        monkeypatch.setattr(options, "can_ask", lambda: can_ask)
    return CliRunner().invoke(cli, ["fetch", "-p", str(root), *args], input=answer)


def _ids(root, name):
    """The `model version` rows of a list, or None without the file."""
    path = root / "input" / f"{name}.csv"
    if not path.exists():
        return None
    lines = path.read_text().split()
    assert lines[0] == "model,version"
    return [line.replace(",", " ") for line in lines[1:]]


def test_fetch_ids_sorts_them_by_role_and_writes_the_lists(tmp_path, fetched):
    root = _inputs(
        tmp_path, [QED, ["ersilia", "eos4zfy", "v1.0.0", "0.5", "higher", "primary"]]
    )

    result = _fetch_ids(root, "eos9p57,eos1bgt", "eos4zfy", "rdkit", "-f")

    assert result.exit_code == 0, result.output
    # a bare id is live
    assert _ids(root, "generators") == ["eos9p57 live"]
    assert _ids(root, "similarity") == ["eos1bgt live"]
    assert fetched["plan"].model_ids == ("eos9p57", "eos1bgt", "eos4zfy")
    # the annotator at the version annotators.csv gives it, and a note says so
    assert fetched["plan"].versions["eos4zfy"] == "v1.0.0"
    assert "lists eos4zfy at v1.0.0, not live" in result.output
    assert fetched["plan"].rdkit
    assert (
        "eos9p57 live, eos1bgt live, eos4zfy v1.0.0, and the rdkit environment"
        in result.output
    )


def test_fetch_id_at_tag_pins_that_release(tmp_path, fetched):
    root = _inputs(
        tmp_path, [QED, ["ersilia", "eos4zfy", "v1.0.0", "0.5", "higher", "primary"]]
    )

    result = _fetch_ids(
        root, "eos9p57@v1.0.0", "eos1bgt@live", "eos4zfy@v2.0.1", "rdkit", "-f"
    )

    assert result.exit_code == 0, result.output
    assert _ids(root, "generators") == ["eos9p57 v1.0.0"]
    assert _ids(root, "similarity") == ["eos1bgt live"]
    # annotators.csv stays the user's: its version is fetched, and the difference said
    assert fetched["plan"].versions["eos4zfy"] == "v1.0.0"
    assert "lists eos4zfy at v1.0.0, not v2.0.1" in result.output


def test_an_empty_tag_after_at_is_refused(tmp_path, fetched):
    result = _fetch_ids(_inputs(tmp_path), "eos9p57@", "rdkit", "-f")

    assert result.exit_code == 1
    assert "a release tag (e.g. eos9p57@v1.0.0) or live" in result.output


def test_a_role_without_ids_has_its_list_removed(tmp_path, fetched):
    root = _inputs(tmp_path)

    result = _fetch_ids(root, "eos9p57", "rdkit", "-f")

    assert result.exit_code == 0, result.output
    assert "similarity.csv: eos8nyg v1.0.0 -> none (file removed)" in result.output
    assert _ids(root, "similarity") is None  # a missing list is none
    _fetch(tmp_path / "models", tmp_path / "envs", "eos9p57")
    inputs = config.load_inputs(str(root / "chemsampler.json"))
    assert inputs.similarity == [] and inputs.generators is not None


@pytest.mark.parametrize(
    "ids, message",
    [
        (["eos9p57", "eos0xxx", "rdkit"], "eos0xxx is not a generator, similarity"),
        (["eos9p57", "eos8nyg"], "no annotator"),
        (["eos4zfy", "rdkit"], "no generator nor similarity model"),
        (["eos9p57", "eos1d7r", "rdkit"], "not meant to be used in chemsampler"),
    ],
)
def test_fetch_ids_that_cannot_make_a_run_are_refused(tmp_path, fetched, ids, message):
    root = _inputs(tmp_path)
    before = (root / "input" / "generators.csv").read_text()

    result = _fetch_ids(root, *ids, "-f")

    assert result.exit_code == 1
    assert message in result.output
    assert "plan" not in fetched
    assert (root / "input" / "generators.csv").read_text() == before


def test_fetch_asks_before_replacing_a_list(tmp_path, fetched, monkeypatch):
    root = _inputs(tmp_path)  # generators: eos9taz, similarity: eos8nyg

    no = _fetch_ids(
        root,
        "eos9p57",
        "eos8nyg@v1.0.0",  # as similarity.csv has it
        "rdkit",
        can_ask=True,
        answer="n\n",
        monkeypatch=monkeypatch,
    )

    assert no.exit_code == 0, no.output
    assert "generators.csv: eos9taz v1.0.0 -> eos9p57 live" in no.output
    assert "similarity.csv" not in no.output  # unchanged: not asked about
    assert "nothing fetched" in no.output and "plan" not in fetched
    assert _ids(root, "generators") == ["eos9taz v1.0.0"]
    yes = _fetch_ids(
        root,
        "eos9p57",
        "eos8nyg@v1.0.0",  # as similarity.csv has it
        "rdkit",
        can_ask=True,
        answer="y\n",
        monkeypatch=monkeypatch,
    )
    assert yes.exit_code == 0, yes.output
    assert _ids(root, "generators") == ["eos9p57 live"]


def test_fetch_cannot_ask_without_a_terminal(tmp_path, fetched, monkeypatch):
    root = _inputs(tmp_path)

    result = _fetch_ids(root, "eos9p57", "rdkit", monkeypatch=monkeypatch)

    assert result.exit_code == 1
    assert "run it with -f" in result.output
    assert _ids(root, "generators") == ["eos9taz v1.0.0"] and "plan" not in fetched


def test_fetch_ids_say_where_annotators_csv_disagrees(tmp_path, fetched):
    root = _inputs(
        tmp_path, [QED, ["ersilia", "eos4zfy", "v1.0.0", "0.5", "higher", "primary"]]
    )

    result = _fetch_ids(root, "eos9p57", "eos9ei3", "-f")

    assert result.exit_code == 0, result.output
    assert (
        "eos9ei3 is not in input/annotators.csv: add a row with its version (live)"
        in result.output
    )
    assert "also lists eos4zfy, which you did not name" in result.output
    assert fetched["plan"].model_ids == ("eos9p57", "eos9ei3")
    assert fetched["plan"].versions["eos9ei3"] == "live"


def test_fetch_does_not_touch_the_store_when_the_inputs_are_wrong(
    tmp_path, store_root, fetched
):
    root = _inputs(tmp_path)
    (root / "input" / "similarity.csv").unlink()
    (root / "input" / "generators.csv").unlink()

    result = CliRunner().invoke(cli, ["fetch", "-p", str(root), "-f"])

    assert result.exit_code == 1
    assert "plan" not in fetched and not store_root.exists()


# --- the project a command works on -------------------------------------------------

#: What each command needs besides the project to succeed here (fakes, rdkit).
NEEDS_A_PROJECT = {
    "fetch": ["-f"],
    "run": ["--seed-smiles", "CCO", "--n-rounds", "1"],
    "annotate": ["--smiles", "CCO"],
}


@pytest.mark.parametrize("command", NEEDS_A_PROJECT)
def test_outside_a_project_a_command_says_to_create_one(
    tmp_path, monkeypatch, fetched, command
):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(cli, [command, *NEEDS_A_PROJECT[command]])

    assert result.exit_code == 1
    assert "not inside a project" in result.output
    assert "chemsampler init" in result.output
    assert "plan" not in fetched


@pytest.mark.parametrize("inside", ["", "input"], ids=["root", "subfolder"])
@pytest.mark.parametrize("command", NEEDS_A_PROJECT)
def test_path_names_the_project_from_anywhere(
    tmp_path, monkeypatch, fetched, command, inside
):
    root = _inputs(tmp_path)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    result = CliRunner().invoke(
        cli, [command, "-p", str(root / inside), *NEEDS_A_PROJECT[command]]
    )

    assert result.exit_code == 0, result.output
    assert list(elsewhere.iterdir()) == []  # results go to the project, not here
    if command == "run":
        assert (_seed_dir(root) / "summary.csv").exists()


@pytest.mark.parametrize("command", [*NEEDS_A_PROJECT, "catalog"])
def test_a_path_in_no_project_is_refused(tmp_path, command):
    plain = tmp_path / "plain"
    plain.mkdir()

    result = CliRunner().invoke(
        cli, [command, "-p", str(plain), *NEEDS_A_PROJECT.get(command, [])]
    )

    assert result.exit_code == 1
    assert f"{plain} is not in a project" in result.output


#: Flags that are gone: where things are is set once, by `init`, and every Hub
#: model runs from its clone (no `--backend`).
REMOVED = {
    "fetch": ["--config", "--path-to-models", "--path-to-model-repositories"]
    + ["--path-to-envs", "--similarity", "--generators", "--annotators", "--rebuild"],
    "run": ["--config", "--path-to-models", "--path-to-model-repositories"]
    + ["--path-to-envs", "--similarity", "--generators", "--annotators"]
    + ["--output-dir", "--backend"],
    "annotate": ["--config", "--path-to-models", "--path-to-model-repositories"]
    + ["--path-to-envs", "--annotators", "--backend"],
    "catalog": ["--config", "--path-to-model-repositories", "--path-to-envs"],
}


@pytest.mark.parametrize(
    "command, flag", [(c, f) for c, flags in REMOVED.items() for f in flags]
)
def test_the_location_flags_are_gone(tmp_path, command, flag):
    result = CliRunner().invoke(cli, [command, flag, str(tmp_path)])

    assert result.exit_code == 2
    assert "No such option" in result.output


# --- the top-level help ---------------------------------------------------------------


def test_every_command_is_described_in_full_in_the_help():
    result = CliRunner().invoke(cli, ["--help"], terminal_width=80)

    assert result.exit_code == 0, result.output
    for name, command in cli.commands.items():
        assert command.short_help, f"{name} has no short description"
        # click cuts a description with "..." when it is too long for the line
        assert command.short_help in result.output, name
    # the usage line ends in `[ARGS]...` by design: look only at the command list
    command_list = result.output.split("Commands:")[1].split("Typical use:")[0]
    assert "..." not in command_list


def test_the_help_lists_commands_in_the_order_a_user_needs_them():
    result = CliRunner().invoke(cli, ["--help"], terminal_width=80)

    positions = [
        result.output.index(f"  {name} ")
        for name in ("init", "fetch", "run", "annotate", "catalog", "suggest")
    ]

    assert positions == sorted(positions)


def test_the_help_ends_with_a_typical_use():
    result = CliRunner().invoke(cli, ["--help"], terminal_width=80)

    for step in ("chemsampler init -p", "chemsampler fetch", "chemsampler run"):
        assert step in result.output
