import sys

import pandas as pd
import pytest

from chemsampler import execution
from chemsampler.hub.client import HubModel
from chemsampler.models.annotator import (
    MolecularWeightAnnotator,
    QEDAnnotator,
    TanimotoAnnotator,
)
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import hill_climb

from .conftest import make_rdkit_env


class RecordingExecutor(execution.LocalExecutor):
    """Runs everything here, as the default does, and records how it was asked."""

    def __init__(self):
        self.label = None
        self.batches = []  # (stage label, number of calls made side by side)
        self.calls = []

    def map(self, function, items):
        items = list(items)
        self.batches.append((self.label, len(items)))
        return super().map(function, items)

    def run_model(self, model, smiles_list):
        self.calls.append(("model", model.model_id, list(smiles_list)))
        return execution.Ran(pd.DataFrame({"smi_0": ["CCN"]}), elapsed=1.5)

    def run_rdkit(self, annotator, smiles_list):
        self.calls.append(("rdkit", annotator.name, list(smiles_list)))
        return super().run_rdkit(annotator, smiles_list)

    def stage(self, label):
        self.label = label
        return super().stage(label)


class Neighbours:
    model_id = "sim"

    def search_by_model(self, seed_smiles):
        return {"sim": ["CCC"]}


class Generator:
    model_id = "gen"

    def generate_by_model(self, seed_smiles):
        return {"gen": ["CCN"]}


class Scores:
    def __init__(self, scores):
        self.scores = scores

    def score(self, smiles_list):
        return {s: self.scores[s] for s in smiles_list if s in self.scores}


def test_each_step_of_a_round_reaches_the_executor_as_one_batch():
    scores = {"CCO": 0.0, "CCN": 1.0, "CCC": 2.0}
    annotators = [
        AnnotatorSpec("a", Scores(scores), cutoff=float("-inf"), direction="higher"),
        AnnotatorSpec("b", Scores(scores), cutoff=float("-inf"), direction="higher"),
    ]
    recorder = RecordingExecutor()

    with execution.use(recorder):
        hill_climb(
            Generator(),
            annotators,
            mode="weighted",
            similarity=Neighbours(),
            seed_smiles="CCO",
            n_rounds=1,
        )

    assert recorder.batches == [
        ("Round 0 · annotation", 2),  # both annotators score the seed together
        ("Round 1 · sampling", 2),  # the generators and the similarity search
        ("Round 1 · annotation", 2),
    ]


def test_hub_models_and_rdkit_annotators_go_through_the_current_executor():
    recorder = RecordingExecutor()
    model = HubModel("eos9taz", models_dir="models")

    with execution.use(recorder):
        table = model.run(["CCO"])
        scores = QEDAnnotator().score(["CCO", "not_a_smiles"])

    assert table["smi_0"].tolist() == ["CCN"] and model.elapsed == 1.5
    assert recorder.calls == [
        ("model", "eos9taz", ["CCO"]),
        ("rdkit", "qed", ["CCO", "not_a_smiles"]),
    ]
    assert list(scores) == ["CCO"]
    assert execution.current() is not recorder  # back to the default afterwards


def test_rdkit_annotators_run_in_the_rdkit_environment_with_the_same_values(tmp_path):
    env = make_rdkit_env(tmp_path / "envs")
    calls = tmp_path / "calls"
    # the env's python notes each call (its folder, its arguments), then runs
    (env / "bin" / "python").write_text(
        f'#!/bin/bash\necho "$PWD $*" >> {calls}\nexec {sys.executable} "$@"\n'
    )
    smiles = ["CCO", "c1ccccc1O", "not_a_smiles"]
    envs = str(tmp_path / "envs")

    for cls, args in [
        (MolecularWeightAnnotator, ()),
        (QEDAnnotator, ()),
        (TanimotoAnnotator, ("CCN",)),
    ]:
        in_env = cls(*args, envs_dir=envs) if args else cls(envs_dir=envs)
        here = cls(*args) if args else cls()
        assert in_env.score(smiles) == here.score(smiles)  # every digit
        assert "not_a_smiles" not in in_env.score(smiles)

    lines = calls.read_text().splitlines()
    assert [line.split()[2] for line in lines] == [
        "mw",
        "mw",
        "qed",
        "qed",
        "tanimoto",
        "tanimoto",
    ]
    assert all(
        "_rdkit_scores.py" in line and str(tmp_path) not in line.split()[0]
        for line in lines
    )  # a temporary folder of its own
    assert lines[-1].endswith("input.csv output.csv CCN")  # the seed


def test_an_rdkit_environment_that_fails_says_so(tmp_path):
    env = make_rdkit_env(tmp_path / "envs")
    (env / "bin" / "python").write_text(
        "#!/bin/bash\necho 'no rdkit here' >&2\nexit 3\n"
    )

    with pytest.raises(
        RuntimeError, match=r"qed: the rdkit environment failed \(exit 3\)"
    ):
        QEDAnnotator(envs_dir=str(tmp_path / "envs")).score(["CCO"])


def test_an_unbuilt_rdkit_environment_says_to_fetch_it(tmp_path):
    with pytest.raises(RuntimeError, match="run `chemsampler fetch rdkit`"):
        QEDAnnotator(envs_dir=str(tmp_path / "envs")).score(["CCO"])
