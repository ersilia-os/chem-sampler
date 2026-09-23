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
    monkeypatch.setattr("chemsampler.hub.client.ErsiliaModel", FakeErsiliaModel)


def test_ersilia_backend_is_default():
    model = HubModel("fake-model")

    assert model.backend == "ersilia"


def test_rejects_unknown_backend():
    with pytest.raises(ValueError, match="backend"):
        HubModel("fake-model", backend="bogus")


def test_ersilia_backend_closes_even_if_serve_fails():
    model = HubModel("fake-model")
    model.model.serve_error = RuntimeError("port in use")

    with pytest.raises(RuntimeError, match="port in use"):
        model.run(["C"])

    assert model.model.closed is True


def test_run_sh_happy_path(tmp_path, monkeypatch):
    bundle = tmp_path
    (bundle / "model" / "framework").mkdir(parents=True)
    (bundle / "service_class.txt").write_text("conda")

    calls = []

    def fake_run(argv, **kwargs):
        calls.append(argv)
        pd.DataFrame({"smiles": ["C"], "score": [1.0]}).to_csv(argv[4], index=False)
        return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

    monkeypatch.setattr("chemsampler.hub.client.subprocess.run", fake_run)

    model = HubModel("fake-model", backend="run_sh")
    model.model.paths = {"repository": str(bundle)}

    df = model.run(["C"])

    assert len(calls) == 1
    argv = calls[0]
    assert argv[0] == "bash"
    assert argv[1] == str(bundle / "model" / "framework" / "run.sh")
    assert argv[2] == str(bundle / "model" / "framework")
    assert argv[5] == str(bundle / "app")
    assert df["score"].tolist() == [1.0]


def test_run_sh_rejects_non_conda_service_class(tmp_path):
    (tmp_path / "service_class.txt").write_text("docker")

    model = HubModel("fake-model", backend="run_sh")
    model.model.paths = {"repository": str(tmp_path)}

    with pytest.raises(RuntimeError, match="conda-packed"):
        model.run(["C"])


def test_run_sh_raises_when_service_class_missing(tmp_path):
    model = HubModel("fake-model", backend="run_sh")
    model.model.paths = {"repository": str(tmp_path)}

    with pytest.raises(RuntimeError, match="service_class.txt"):
        model.run(["C"])


def test_run_sh_raises_when_bundle_unresolved():
    model = HubModel("fake-model", backend="run_sh")
    model.model.paths = {"repository": None}

    with pytest.raises(RuntimeError, match="not fetched"):
        model.run(["C"])


def test_run_sh_wraps_subprocess_failure(tmp_path, monkeypatch):
    (tmp_path / "service_class.txt").write_text("conda")

    def fake_run(argv, **kwargs):
        raise subprocess.CalledProcessError(1, argv, output="", stderr="boom")

    monkeypatch.setattr("chemsampler.hub.client.subprocess.run", fake_run)

    model = HubModel("fake-model", backend="run_sh")
    model.model.paths = {"repository": str(tmp_path)}

    with pytest.raises(RuntimeError, match="boom"):
        model.run(["C"])
