import pytest

from chemsampler import store


def test_everything_is_inside_one_folder(store_root):
    resolved = store.resolve()

    assert resolved.root == store_root
    assert resolved.models == store_root / "models"
    assert resolved.envs == store_root / "envs"
    assert resolved.logs == store_root / "logs"
    assert not store_root.exists()  # nothing is created just by asking


@pytest.mark.parametrize("unset", ["delete", "empty"])
def test_the_folder_defaults_to_dot_chemsampler_in_home(tmp_path, monkeypatch, unset):
    monkeypatch.setenv("HOME", str(tmp_path / "me"))
    if unset == "delete":
        monkeypatch.delenv(store.HOME_VARIABLE)
    else:
        monkeypatch.setenv(store.HOME_VARIABLE, "")

    assert store.resolve().root == tmp_path / "me" / ".chemsampler"


def test_a_home_shortcut_in_the_variable_is_expanded(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "me"))
    monkeypatch.setenv(store.HOME_VARIABLE, "~/shared/.chemsampler")

    assert store.home() == tmp_path / "me" / "shared" / ".chemsampler"


def test_an_old_config_file_is_refused_not_ignored(store_root):
    # it used to move the store: ignoring it would make every model look unfetched
    store_root.mkdir()
    (store_root / "config.json").write_text('{"store": "/somewhere/else"}')

    with pytest.raises(ValueError, match="no longer read.*CHEMSAMPLER_HOME"):
        store.resolve()


def test_what_fetch_downloaded_is_recorded_beside_the_clones(tmp_path):
    # never in the store root: two projects may keep their clones apart, and the
    # record describes one set of clones
    assert store.state_path(tmp_path / "models") == (
        tmp_path / "models" / ".chemsampler-fetch.json"
    )
