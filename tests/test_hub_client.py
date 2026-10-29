import subprocess

import pandas as pd
import pytest

from chemsampler.hub.client import HubModel


class FakeErsiliaModel:
    """Stands in for ersilia.ErsiliaModel, matching its real construction shape."""

    def __init__(self, model):
        self.model_id = model
        self.paths = {"repository": None}
        self.closed = False
        self.serve_error = None

    def serve(self):
        if self.serve_error is not None:
            raise self.serve_error

    def run(self, input, output):
        pd.DataFrame({"smiles": ["C"], "score": [1.0]}).to_csv(output, index=False)

    def close(self):
        self.closed = True


@pytest.fixture(autouse=True)
def fake_ersilia_model(monkeypatch):
    monkeypatch.setattr("ersilia.ErsiliaModel", FakeErsiliaModel)


def test_run_sh_backend_is_default():
    model = HubModel("fake-model", models_dir="/models")

    assert model.backend == "run_sh"


def test_run_sh_requires_the_models_dir():
    with pytest.raises(ValueError, match="path-to-model-repositories"):
        HubModel("fake-model")


def test_ersilia_backend_does_not_need_the_models_dir():
    assert HubModel("fake-model", backend="ersilia").models_dir is None


def test_rejects_unknown_backend():
    with pytest.raises(ValueError, match="backend"):
        HubModel("fake-model", backend="bogus")


def test_ersilia_backend_closes_even_if_serve_fails():
    model = HubModel("fake-model", backend="ersilia")
    model.model.serve_error = RuntimeError("port in use")

    with pytest.raises(RuntimeError, match="port in use"):
        model.run(["C"])

    assert model.model.closed is True


def test_run_sh_happy_path(tmp_path, monkeypatch):
    models_dir = tmp_path / "models"
    bundle = models_dir / "fake-model"
    (bundle / "model" / "framework").mkdir(parents=True)
    (bundle / "app").mkdir(parents=True)
    (bundle / "model" / "framework" / "run.sh").write_text("#!/bin/bash\necho ok")

    calls = []

    def fake_run(argv, **kwargs):
        calls.append(argv)
        pd.DataFrame({"smiles": ["C"], "score": [1.0]}).to_csv(argv[4], index=False)
        return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

    monkeypatch.setattr("chemsampler.hub.client.subprocess.run", fake_run)

    model = HubModel("fake-model", backend="run_sh", models_dir=str(models_dir))
    df = model.run(["C"])

    assert len(calls) == 1
    argv = calls[0]
    assert argv[0] == "bash"
    assert argv[1] == str(bundle / "model" / "framework" / "run.sh")
    assert argv[2] == str(bundle / "model" / "framework")
    assert argv[5] == str(bundle / "app")
    assert df["score"].tolist() == [1.0]


def test_run_sh_raises_when_run_sh_missing(tmp_path):
    models_dir = tmp_path / "models"
    bundle = models_dir / "fake-model"
    (bundle / "model" / "framework").mkdir(parents=True)

    model = HubModel("fake-model", backend="run_sh", models_dir=str(models_dir))

    with pytest.raises(RuntimeError, match="run.sh not found"):
        model.run(["C"])


def test_run_sh_raises_when_bundle_missing():
    model = HubModel("fake-model", backend="run_sh", models_dir="/nonexistent")

    with pytest.raises(RuntimeError, match="not found at"):
        model.run(["C"])


def test_run_sh_wraps_subprocess_failure(tmp_path, monkeypatch):
    models_dir = tmp_path / "models"
    bundle = models_dir / "fake-model"
    (bundle / "model" / "framework").mkdir(parents=True)
    (bundle / "model" / "framework" / "run.sh").write_text("#!/bin/bash")

    def fake_run(argv, **kwargs):
        raise subprocess.CalledProcessError(1, argv, output="", stderr="boom")

    monkeypatch.setattr("chemsampler.hub.client.subprocess.run", fake_run)

    model = HubModel("fake-model", backend="run_sh", models_dir=str(models_dir))

    with pytest.raises(RuntimeError, match="boom"):
        model.run(["C"])
