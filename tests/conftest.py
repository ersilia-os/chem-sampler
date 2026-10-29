import pytest

from chemsampler import store


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
