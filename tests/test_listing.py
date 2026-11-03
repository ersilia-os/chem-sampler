import subprocess

import pytest

from chemsampler.hub.catalog import Catalog, CatalogEntry
from chemsampler.hub.envs import env_path, write_marker
from chemsampler.hub.listing import (
    env_state,
    list_hub,
    list_local,
    present_models,
    version_at,
)

ENTRIES = (
    CatalogEntry("generators", "eos9p57", "CReM Growth", "Sampling", "Generation"),
    CatalogEntry("generators", "eos9taz", "MoLeR", "Sampling", "Generation"),
    CatalogEntry(
        "similarity", "eos8nyg", "Reference library", "Sampling", "Similarity search"
    ),
    CatalogEntry(
        "similarity", "eos2a9n", "ChEMBL search", "Sampling", "Similarity search"
    ),
    CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", ""),
)
CATALOG = Catalog("2026-10-29", "0.1.62", ENTRIES)


@pytest.fixture
def folders(tmp_path):
    return tmp_path / "models", tmp_path / "envs"


def _clone(models, model_id):
    (models / model_id).mkdir(parents=True)


def _env(envs, model_id, complete=True):
    env_path(envs, model_id).mkdir(parents=True)
    if complete:
        write_marker(env_path(envs, model_id), model=model_id)


# --- what is on disk ----------------------------------------------------------------


def test_clones_are_the_folders_and_nothing_else(folders):
    models, _ = folders
    _clone(models, "eos9taz")
    (models / ".chemsampler-fetch.json").write_text("{}")  # fetch's own record
    (models / "notes.txt").write_text("x")
    (models / ".cache").mkdir()

    assert present_models(models) == {"eos9taz"}


def test_a_missing_models_folder_is_simply_empty(tmp_path):
    assert present_models(tmp_path / "nowhere") == set()
    (tmp_path / "file").write_text("x")
    assert present_models(tmp_path / "file") == set()


def test_an_environment_is_ok_unfinished_or_missing(folders):
    _, envs = folders
    _env(envs, "eos9taz")
    _env(envs, "eos8nyg", complete=False)

    assert env_state(envs, "eos9taz") == "ok"
    assert env_state(envs, "eos8nyg") == "unfinished"
    assert env_state(envs, "eos4zfy") == "missing"


def test_listing_never_reads_inside_an_environment(folders, monkeypatch):
    # `pip list` and `du` took tens of seconds over NFS: only the marker is read
    models, envs = folders
    _clone(models, "eos9taz")
    _env(envs, "eos9taz")
    opened = []
    real_read = type(models).read_text

    def spy(self, *args, **kwargs):
        opened.append(self.name)
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(type(models), "read_text", spy)

    list_local(models, envs, catalog=CATALOG)

    assert set(opened) <= {".chemsampler-env.json"}  # the marker, nothing else


# --- the local listing --------------------------------------------------------------


def test_each_model_gets_its_role_title_and_environment_state(folders):
    models, envs = folders
    for model_id in ("eos9p57", "eos8nyg", "eos4zfy"):
        _clone(models, model_id)
    _env(envs, "eos9p57")
    _env(envs, "eos8nyg", complete=False)

    listed = {m.model_id: m for m in list_local(models, envs, catalog=CATALOG)}

    assert (listed["eos9p57"].role, listed["eos9p57"].title) == (
        "generators",
        "CReM Growth",
    )
    assert listed["eos9p57"].env == "ok"
    assert listed["eos8nyg"].env == "unfinished"
    assert listed["eos4zfy"].env == "missing" and listed["eos4zfy"].cloned


def test_models_are_ordered_by_role_then_id(folders):
    models, envs = folders
    for model_id in ("eos4zfy", "eos9taz", "eos8nyg", "eos9p57"):
        _clone(models, model_id)

    ids = [m.model_id for m in list_local(models, envs, catalog=CATALOG)]

    assert ids == ["eos9p57", "eos9taz", "eos8nyg", "eos4zfy"]


def test_a_model_the_project_needs_but_the_store_lacks_is_listed(folders):
    models, envs = folders
    _clone(models, "eos9taz")

    listed = {
        m.model_id: m
        for m in list_local(models, envs, used=["eos9taz", "eos8nyg"], catalog=CATALOG)
    }

    assert listed["eos8nyg"].used and not listed["eos8nyg"].cloned
    assert listed["eos8nyg"].env == "missing"
    assert listed["eos8nyg"].problem == "not cloned"
    assert listed["eos9taz"].used and listed["eos9taz"].cloned


def test_a_model_the_project_does_not_need_is_not_marked(folders):
    models, envs = folders
    _clone(models, "eos9taz")

    [model] = list_local(models, envs, used=[], catalog=CATALOG)

    assert not model.used and model.problem is None


def test_an_unknown_folder_is_listed_not_hidden(folders):
    models, envs = folders
    _clone(models, "my-experiment")

    [model] = list_local(models, envs, catalog=CATALOG)

    assert (model.role, model.title) == ("unknown", "")


def test_without_a_catalog_every_role_is_unknown(folders):
    models, envs = folders
    _clone(models, "eos9taz")

    [model] = list_local(models, envs, catalog=None)

    assert model.role == "unknown"


def test_no_environments_folder_means_no_environment_to_report(folders):
    # models fetched by `ersilia fetch` bring their own interpreter
    models, _ = folders
    _clone(models, "eos9taz")

    [model] = list_local(models, None, catalog=CATALOG)

    assert model.env == "-"


def test_the_local_listing_can_be_filtered_by_role(folders):
    models, envs = folders
    for model_id in ("eos9taz", "eos8nyg", "eos4zfy"):
        _clone(models, model_id)

    assert [
        m.model_id for m in list_local(models, envs, catalog=CATALOG, role="similarity")
    ] == ["eos8nyg"]


def test_an_empty_store_lists_nothing(folders):
    models, envs = folders

    assert list_local(models, envs, catalog=CATALOG) == []


# --- the Hub listing ----------------------------------------------------------------


def test_the_hub_listing_marks_what_this_machine_has(folders):
    models, _ = folders
    _clone(models, "eos9taz")

    listed = {e.model_id: e for e in list_hub(CATALOG, models)}

    assert listed["eos9taz"].in_store and not listed["eos9p57"].in_store
    assert len(listed) == 5


def test_only_the_accepted_similarity_models_are_supported(folders):
    models, _ = folders

    listed = {e.model_id: e for e in list_hub(CATALOG, models, role="similarity")}

    assert listed["eos8nyg"].supported and not listed["eos2a9n"].supported
    assert set(listed) == {"eos8nyg", "eos2a9n"}


def test_every_generator_and_annotator_counts_as_supported(folders):
    models, _ = folders

    assert all(
        e.supported
        for e in list_hub(CATALOG, models)
        if e.role in ("generators", "annotators")
    )


# --- versions and stale environments ------------------------------------------------


def _git(folder, *args):
    out = subprocess.run(
        ["git", "-C", str(folder), "-c", "user.name=t", "-c", "user.email=t@t", *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout.strip()


def _git_clone(models, model_id):
    """A clone on its main branch with a v1.0.0 release and one commit after it."""
    folder = models / model_id
    folder.mkdir(parents=True)
    _git(folder, "init", "-q", "-b", "main")
    (folder / "a").write_text("1")
    _git(folder, "add", "a")
    _git(folder, "commit", "-q", "-m", "first")
    _git(folder, "tag", "v1.0.0")
    (folder / "a").write_text("2")
    _git(folder, "commit", "-q", "-am", "more work")
    return folder


def test_a_clone_on_the_wanted_version_shows_it(folders):
    models, envs = folders
    _git_clone(models, "eos9taz")

    assert version_at(models / "eos9taz", "live") == "live"
    [model] = list_local(
        models, envs, used=["eos9taz"], versions={"eos9taz": "live"}, catalog=CATALOG
    )
    assert (model.version, model.at) == ("live", "live")
    assert "the inputs ask for" not in (model.problem or "")


def test_a_clone_at_another_version_than_the_inputs_ask_for_says_so(folders):
    models, envs = folders
    folder = _git_clone(models, "eos9taz")
    _git(folder, "checkout", "-q", "v1.0.0")

    [model] = list_local(
        models, envs, used=["eos9taz"], versions={"eos9taz": "live"}, catalog=CATALOG
    )

    assert model.at == "v1.0.0"
    assert "cloned at v1.0.0, the inputs ask for live" in model.problem


def test_a_clone_on_no_branch_nor_tag_shows_its_commit(folders):
    models, _ = folders
    folder = _git_clone(models, "eos9taz")
    _git(folder, "checkout", "-q", "--detach", "HEAD~1")
    _git(folder, "tag", "-d", "v1.0.0")

    assert version_at(folder) == _git(folder, "rev-parse", "--short", "HEAD")


def test_a_folder_that_is_no_git_clone_has_no_version_to_compare(folders):
    models, envs = folders
    _clone(models, "eos9taz")

    [model] = list_local(
        models, envs, used=["eos9taz"], versions={"eos9taz": "v9.9.9"}, catalog=CATALOG
    )

    assert model.at is None and model.problem == "environment missing"


def test_an_environment_built_from_another_commit_is_stale(folders):
    models, envs = folders
    folder = _git_clone(models, "eos9taz")
    built = _git(folder, "rev-parse", "--short", "HEAD")
    env_path(envs, "eos9taz").mkdir(parents=True)
    write_marker(env_path(envs, "eos9taz"), model="eos9taz", commit=built)
    assert env_state(envs, "eos9taz", models) == "ok"

    _git(folder, "checkout", "-q", "v1.0.0")

    assert env_state(envs, "eos9taz", models) == "stale"
    [model] = list_local(
        models, envs, used=["eos9taz"], versions={"eos9taz": "v1.0.0"}, catalog=CATALOG
    )
    head = _git(folder, "rev-parse", "--short", "HEAD")
    assert model.env == "stale"
    assert model.problem == f"environment built from {built}, the clone is at {head}"


def test_an_environment_that_records_no_commit_is_never_stale(folders):
    models, envs = folders
    _git_clone(models, "eos9taz")
    _env(envs, "eos9taz")  # a marker without a commit

    assert env_state(envs, "eos9taz", models) == "ok"
