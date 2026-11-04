import json
import shutil
import subprocess

import pytest
from click.testing import CliRunner

from chemsampler import config
from chemsampler.cli import cli, options
from chemsampler.hub import remove as hub_remove
from chemsampler.store import STATE_FILE

from .conftest import make_project, make_rdkit_env


@pytest.fixture(autouse=True)
def _no_hub_check(monkeypatch):
    monkeypatch.setattr(config, "validate_models", lambda **kwargs: None)


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A project whose template uses eos9taz, eos8nyg and eos4zfy (all live)."""
    make_project(tmp_path / "proj")
    monkeypatch.chdir(tmp_path / "proj")
    return tmp_path / "proj"


def _remove(*args, ask=False, answer=None, monkeypatch=None):
    if monkeypatch is not None:
        monkeypatch.setattr(options, "can_ask", lambda: ask)
    return CliRunner().invoke(cli, ["remove", *args], input=answer)


def _there(store_root, model_id):
    """(clone, environment) of a model on disk."""
    return (
        (store_root / "models" / model_id).exists(),
        (store_root / "envs" / model_id).exists(),
    )


def _git(folder, *args):
    subprocess.run(
        ["git", "-C", str(folder), "-c", "user.name=t", "-c", "user.email=t@t", *args],
        capture_output=True,
        text=True,
        check=True,
    )


def _git_clone(store_root, model_id):
    """A clone whose HEAD is a release: nothing only this clone has."""
    folder = store_root / "models" / model_id
    folder.mkdir(parents=True)
    _git(folder, "init", "-q", "-b", "main")
    (folder / "a").write_text("1")
    _git(folder, "add", "a")
    _git(folder, "commit", "-q", "-m", "first")
    _git(folder, "tag", "v1.0.0")
    return folder


def _record(store_root, **entries):
    (store_root / "models").mkdir(parents=True, exist_ok=True)
    (store_root / "models" / STATE_FILE).write_text(json.dumps(entries))


def _recorded(store_root):
    return json.loads((store_root / "models" / STATE_FILE).read_text())


# --- what it removes ----------------------------------------------------------------


def test_remove_deletes_the_clone_the_environment_and_the_record(
    project, store_root, fake_fetch
):
    fake_fetch("eos9taz", "eos9p57")
    _record(store_root, eos9taz={"checkpoints_commit": "x"}, eos9p57={"c": 1})
    csvs = {p.name: p.read_bytes() for p in (project / "input").glob("*.csv")}

    result = _remove("eos9p57", "-f")

    assert result.exit_code == 0, result.output
    assert _there(store_root, "eos9p57") == (False, False)
    assert _there(store_root, "eos9taz") == (True, True)  # the others stay
    assert _recorded(store_root) == {"eos9taz": {"checkpoints_commit": "x"}}
    assert "Removed 1 model(s): eos9p57" in result.output
    assert csvs == {p.name: p.read_bytes() for p in (project / "input").glob("*.csv")}


def test_ids_can_be_separated_by_commas(project, store_root, fake_fetch):
    fake_fetch("eos9p57", "eos1bgt")

    result = _remove("eos9p57,eos1bgt", "-f")

    assert result.exit_code == 0, result.output
    assert not any(_there(store_root, m)[0] for m in ("eos9p57", "eos1bgt"))


def test_an_unfinished_environment_a_lone_environment_and_a_lone_clone_all_go(
    project, store_root, fake_fetch
):
    fake_fetch("eos9p57", complete=False)  # clone + unfinished env
    (store_root / "envs" / "eos1bgt").mkdir(parents=True)  # env without a clone
    (store_root / "models" / "eos8fma").mkdir(parents=True)  # clone without an env

    result = _remove("eos9p57", "eos1bgt", "eos8fma", "-f")

    assert result.exit_code == 0, result.output
    for model_id in ("eos9p57", "eos1bgt", "eos8fma"):
        assert _there(store_root, model_id) == (False, False)


def test_a_model_that_is_not_in_the_store_is_skipped(project, store_root, fake_fetch):
    fake_fetch("eos9p57")

    result = _remove("eos9p57", "eos0xxx", "-f")

    assert result.exit_code == 0, result.output
    assert "eos0xxx: not in the store" in result.output
    assert _there(store_root, "eos9p57") == (False, False)


def test_nothing_in_the_store_to_remove_is_an_error(project, store_root):
    result = _remove("eos0xxx", "-f")

    assert result.exit_code == 1
    assert "eos0xxx: not in the store" in result.output
    assert "Nothing removed" in result.output


def test_remove_rdkit_deletes_that_environment_only(project, store_root, fake_fetch):
    make_rdkit_env(store_root / "envs")
    fake_fetch("eos9taz")

    result = _remove("rdkit", "-f")

    assert result.exit_code == 0, result.output
    assert not (store_root / "envs" / "rdkit").exists()
    assert _there(store_root, "eos9taz") == (True, True)


def test_remove_rdkit_never_touches_a_clone_folder_of_that_name(project, store_root):
    make_rdkit_env(store_root / "envs")
    (store_root / "models" / "rdkit").mkdir(parents=True)  # the user's own folder

    result = _remove("rdkit", "-f")

    assert result.exit_code == 0, result.output
    assert not (store_root / "envs" / "rdkit").exists()
    assert (store_root / "models" / "rdkit").exists()


@pytest.mark.parametrize("name", ["mw", "qed", "tanimoto"])
def test_the_rdkit_annotators_point_to_remove_rdkit(project, store_root, name):
    result = _remove(name, "-f")

    assert result.exit_code == 1
    assert "chemsampler remove rdkit" in result.output


# --- asking -------------------------------------------------------------------------


def test_it_lists_what_goes_and_asks_before_removing(
    project, store_root, fake_fetch, monkeypatch
):
    fake_fetch("eos9taz", "eos9p57")

    refused = _remove(
        "eos9taz", "eos9p57", ask=True, answer="n\n", monkeypatch=monkeypatch
    )

    assert refused.exit_code == 0, refused.output
    assert "eos9p57  clone, environment" in refused.output
    # eos9taz is in the project's generators.csv, eos9p57 is not
    assert "eos9taz is used by this project" in refused.output
    assert "eos9p57 is used" not in refused.output
    assert "Other projects sharing this store" in refused.output
    assert "Nothing removed" in refused.output
    assert _there(store_root, "eos9taz") == (True, True)

    agreed = _remove("eos9p57", ask=True, answer="y\n", monkeypatch=monkeypatch)

    assert agreed.exit_code == 0, agreed.output
    assert _there(store_root, "eos9p57") == (False, False)


def test_without_a_terminal_it_needs_force(
    project, store_root, fake_fetch, monkeypatch
):
    fake_fetch("eos9p57")

    result = _remove("eos9p57", ask=False, monkeypatch=monkeypatch)

    assert result.exit_code == 1
    assert "run it with -f" in result.output
    assert _there(store_root, "eos9p57") == (True, True)


def test_unreadable_inputs_only_drop_the_used_by_this_project_line(
    project, store_root, fake_fetch
):
    fake_fetch("eos9taz")
    (project / "input" / "annotators.csv").write_text("garbage\n")

    result = _remove("eos9taz", "-f")

    assert result.exit_code == 0, result.output
    assert "could not be read" in result.output
    assert _there(store_root, "eos9taz") == (False, False)


# --- local work ---------------------------------------------------------------------


def test_a_clone_with_local_changes_is_kept(project, store_root):
    folder = _git_clone(store_root, "eos9p57")
    (folder / "a").write_text("edited")

    result = _remove("eos9p57")  # nothing else to remove, so nothing to ask

    assert result.exit_code == 1
    assert "eos9p57: local changes, not removed (use -f" in result.output
    assert "Nothing removed" in result.output
    assert (folder / "a").read_text() == "edited"


def test_local_work_blocks_one_model_and_not_the_others(
    project, store_root, fake_fetch, monkeypatch
):
    folder = _git_clone(store_root, "eos9p57")
    (folder / "a").write_text("edited")
    fake_fetch("eos1bgt")

    result = _remove(
        "eos9p57", "eos1bgt", ask=True, answer="y\n", monkeypatch=monkeypatch
    )

    assert result.exit_code == 1  # something asked for was kept
    assert "eos9p57: local changes, not removed (use -f" in result.output
    assert _there(store_root, "eos1bgt") == (False, False)
    assert (store_root / "models" / "eos9p57").exists()


def test_local_commits_also_block_removal(project, store_root):
    folder = _git_clone(store_root, "eos9p57")
    (folder / "a").write_text("2")
    _git(folder, "commit", "-q", "-am", "mine")  # after the release: only here

    result = _remove("eos9p57")

    assert result.exit_code == 1
    assert "eos9p57: local commits, not removed" in result.output
    assert (store_root / "models" / "eos9p57").exists()


def test_force_removes_a_clone_that_holds_local_work(project, store_root):
    folder = _git_clone(store_root, "eos9p57")
    (folder / "a").write_text("edited")

    result = _remove("eos9p57", "-f")

    assert result.exit_code == 0, result.output
    assert not (store_root / "models" / "eos9p57").exists()


def test_a_clean_release_clone_is_removed(project, store_root):
    _git_clone(store_root, "eos9p57")

    assert _remove("eos9p57", "-f").exit_code == 0
    assert not (store_root / "models" / "eos9p57").exists()


def test_a_folder_that_is_no_git_clone_is_removed_without_the_check(
    project, store_root, monkeypatch
):
    (store_root / "models" / "eos9p57").mkdir(parents=True)
    (store_root / "models" / "eos9p57" / "x").write_text("x")

    # no -f, and still not refused: it has no local work to lose
    result = _remove("eos9p57", ask=True, answer="y\n", monkeypatch=monkeypatch)

    assert result.exit_code == 0, result.output
    assert not (store_root / "models" / "eos9p57").exists()


# --- --unused -----------------------------------------------------------------------


def test_unused_removes_what_the_project_does_not_name(project, store_root, fake_fetch):
    fake_fetch("eos9taz", "eos8nyg", "eos4zfy")  # the project's three
    fake_fetch("eos9p57", "eos1bgt")  # spare
    (store_root / "envs" / "eos8fma").mkdir(parents=True)  # an orphan environment
    make_rdkit_env(store_root / "envs")  # no rdkit rows in the template

    result = _remove("--unused", "-f")

    assert result.exit_code == 0, result.output
    for kept in ("eos9taz", "eos8nyg", "eos4zfy"):
        assert _there(store_root, kept) == (True, True)
    for gone in ("eos9p57", "eos1bgt", "eos8fma"):
        assert _there(store_root, gone) == (False, False)
    assert not (store_root / "envs" / "rdkit").exists()


def test_unused_keeps_the_rdkit_environment_when_an_annotator_needs_it(
    project, store_root, fake_fetch
):
    (project / "input" / "annotators.csv").write_text(
        "source,model,version,cutoff,direction,priority\n"
        "ersilia,eos4zfy,live,0.5,higher,primary\n"
        "rdkit,qed,,0.5,higher,guard-rail\n"
    )
    make_rdkit_env(store_root / "envs")
    fake_fetch("eos9taz", "eos8nyg", "eos4zfy")

    result = _remove("--unused", "-f")

    assert result.exit_code == 0, result.output
    assert "Nothing unused" in result.output
    assert (store_root / "envs" / "rdkit").exists()


def test_unused_with_nothing_unused_is_not_an_error(project, store_root, fake_fetch):
    fake_fetch("eos9taz", "eos8nyg", "eos4zfy")

    result = _remove("--unused", "-f")

    assert result.exit_code == 0 and "Nothing unused" in result.output


def test_unused_never_touches_the_record_file_or_hidden_folders(
    project, store_root, fake_fetch
):
    fake_fetch("eos9taz", "eos8nyg", "eos4zfy")
    _record(store_root, eos9taz={"c": 1})
    (store_root / "models" / ".cache").mkdir()

    assert _remove("--unused", "-f").exit_code == 0
    assert (store_root / "models" / STATE_FILE).exists()
    assert (store_root / "models" / ".cache").exists()


def test_unused_needs_readable_inputs(project, store_root, fake_fetch):
    fake_fetch("eos9p57")
    (project / "input" / "annotators.csv").write_text("garbage\n")

    result = _remove("--unused", "-f")

    assert result.exit_code == 1
    assert "could not be read" in result.output
    assert _there(store_root, "eos9p57") == (True, True)


def test_unused_respects_local_work(project, store_root, fake_fetch):
    folder = _git_clone(store_root, "eos9p57")
    (folder / "a").write_text("edited")
    fake_fetch("eos9taz", "eos8nyg", "eos4zfy")

    result = _remove("--unused")

    assert result.exit_code == 1
    assert "eos9p57: local changes, not removed" in result.output
    assert folder.exists()


# --- usage, projects, failures ------------------------------------------------------


@pytest.mark.parametrize("args", [(), ("eos9p57", "--unused")])
def test_it_takes_ids_or_unused_but_not_both_nor_neither(project, args):
    result = _remove(*args)

    assert result.exit_code == 2
    assert "ids, or --unused" in result.output


def test_outside_a_project_it_refuses(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = _remove("eos9p57", "-f")

    assert result.exit_code == 1
    assert "not inside a project" in result.output


def test_a_project_can_be_named_with_p(
    project, store_root, fake_fetch, tmp_path, monkeypatch
):
    fake_fetch("eos9p57")
    monkeypatch.chdir(tmp_path)

    result = _remove("eos9p57", "-f", "-p", str(project))

    assert result.exit_code == 0, result.output
    assert _there(store_root, "eos9p57") == (False, False)


def test_one_failure_does_not_stop_the_others(
    project, store_root, fake_fetch, monkeypatch
):
    fake_fetch("eos9p57", "eos1bgt")
    real = shutil.rmtree

    def flaky(path, *args, **kwargs):
        if "eos9p57" in str(path):
            raise PermissionError("read-only")
        return real(path, *args, **kwargs)

    monkeypatch.setattr(hub_remove.shutil, "rmtree", flaky)

    result = _remove("eos9p57", "eos1bgt", "-f")

    assert result.exit_code == 1
    assert "eos9p57: could not be removed: read-only" in result.output
    assert _there(store_root, "eos1bgt") == (False, False)
    assert "Removed 1 model(s): eos1bgt" in result.output


def test_without_a_record_file_none_is_made(project, store_root, fake_fetch):
    fake_fetch("eos9p57")  # no record file at all

    assert _remove("eos9p57", "-f").exit_code == 0
    assert not (store_root / "models" / STATE_FILE).exists()


def test_remove_is_listed_among_the_commands():
    assert "remove" in CliRunner().invoke(cli, ["--help"]).output
