import pytest

from chemsampler import store
from chemsampler.hub.envs import write_marker


@pytest.fixture(autouse=True)
def store_root(tmp_path, monkeypatch):
    """
    Point the store somewhere harmless for every test.

    Without this a test would read, and `fetch` would write, the real
    `~/.chemsampler` of whoever runs the suite.
    """
    root = tmp_path / "store"
    monkeypatch.setattr(store, "DEFAULT_STORE", root)
    return root


@pytest.fixture
def fake_fetch(store_root):
    """
    Make models look fetched into the test store, as `chemsampler fetch` leaves
    them: a clone folder and an environment, finished unless `complete` is False.
    """

    def fetch(*model_ids, complete=True):
        for model_id in model_ids:
            (store_root / "models" / model_id).mkdir(parents=True)
            env = store_root / "envs" / model_id
            env.mkdir(parents=True)
            if complete:
                write_marker(env, model=model_id)

    return fetch
