import pandas as pd

from chemsampler import execution
from chemsampler.hub.client import HubModel
from chemsampler.models.annotator import QEDAnnotator
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import hill_climb


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

    def run_rdkit(self, name, smiles_list):
        self.calls.append(("rdkit", name, list(smiles_list)))
        return super().run_rdkit(name, smiles_list)

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
