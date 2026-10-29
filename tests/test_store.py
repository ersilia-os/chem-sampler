import json

import pytest

from chemsampler import store


def _write_config(**keys):
    store.config_path().parent.mkdir(parents=True, exist_ok=True)
    store.config_path().write_text(json.dumps(keys))


def test_the_store_defaults_to_one_folder_per_machine(store_root):
    resolved = store.resolve()

    assert resolved.root == store_root
    assert resolved.models == store_root / "models"
    assert resolved.envs == store_root / "envs"
    assert resolved.logs == store_root / "logs"


def test_no_settings_file_is_normal(store_root):
    assert store.read_config() == {}
    assert not store_root.exists()  # nothing is created just by asking


def test_the_whole_store_can_be_moved(store_root, tmp_path):
    _write_config(store=str(tmp_path / "shared"))

    resolved = store.resolve()

    assert resolved.root == tmp_path / "shared"
    assert resolved.models == tmp_path / "shared" / "models"
    assert resolved.envs == tmp_path / "shared" / "envs"


def test_a_single_folder_can_be_moved(store_root, tmp_path):
    _write_config(path_to_envs=str(tmp_path / "fast-disk" / "envs"))

    resolved = store.resolve()

    assert resolved.models == store_root / "models"  # untouched
    assert resolved.envs == tmp_path / "fast-disk" / "envs"


def test_a_relative_override_is_read_inside_the_store(store_root):
    _write_config(path_to_model_repositories="clones")

    assert store.resolve().models == store_root / "clones"


def test_a_home_shortcut_is_expanded(store_root, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "elsewhere"))
    _write_config(store="~/shared-store")

    assert store.resolve().root == tmp_path / "elsewhere" / "shared-store"


@pytest.mark.parametrize(
    "text, message",
    [
        ("[1, 2]", "expected a JSON object"),
        ("{oops", "not valid JSON"),
        ('{"colour": "red"}', "unknown key"),
        ('{"store": 3}', "must be a path"),
    ],
)
def test_broken_settings_are_explained(store_root, text, message):
    store.config_path().parent.mkdir(parents=True, exist_ok=True)
    store.config_path().write_text(text)

    with pytest.raises(ValueError, match=message):
        store.resolve()


def test_the_settings_file_sits_next_to_the_default_store(store_root):
    assert store.config_path() == store_root / "config.json"


def test_what_fetch_downloaded_is_recorded_beside_the_clones(tmp_path):
    # never in the store root: two projects may keep their clones apart, and the
    # record describes one set of clones
    assert store.state_path(tmp_path / "models") == (
        tmp_path / "models" / ".chemsampler-fetch.json"
    )
