import json
import re

import pytest
from click.testing import CliRunner

from chemsampler import config
from chemsampler.cli import cli
from chemsampler.hub import catalog
from chemsampler.hub.catalog import CatalogEntry, write_catalog
from chemsampler.hub.envs import write_marker
from chemsampler.project import create_project

ENTRIES = [
    CatalogEntry("generators", "eos9p57", "CReM Growth", "Sampling", "Generation"),
    CatalogEntry("generators", "eos9taz", "MoLeR", "Sampling", "Generation"),
    CatalogEntry(
        "similarity", "eos8nyg", "Reference library", "Sampling", "Similarity search"
    ),
    CatalogEntry(
        "similarity", "eos2a9n", "ChEMBL search", "Sampling", "Similarity search"
    ),
    CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", ""),
]


@pytest.fixture(autouse=True)
def hub_catalog(tmp_path, monkeypatch):
    """A small catalog file, and no Hub check of the project's models."""
    path = tmp_path / "hub-catalog.csv"
    write_catalog(ENTRIES, path)
    monkeypatch.setattr(catalog, "CATALOG_PATH", path)
    monkeypatch.setattr(config, "validate_models", lambda **kwargs: None)
    return path


def _fetched(store_root, *model_ids, complete=True):
    for model_id in model_ids:
        (store_root / "models" / model_id).mkdir(parents=True)
        env = store_root / "envs" / model_id
        env.mkdir(parents=True)
        if complete:
            write_marker(env, model=model_id)


def _catalog(*args):
    return CliRunner().invoke(cli, ["catalog", *args])


# --- what this machine has -------------------------------------------------------------


def test_the_local_listing_says_where_it_looked_and_what_it_found(tmp_path, store_root):
    _fetched(store_root, "eos9p57", "eos8nyg", "eos4zfy")

    result = _catalog()

    assert result.exit_code == 0, result.output
    assert f"Models in {store_root / 'models'}" in result.output
    assert f"environments in {store_root / 'envs'}" in result.output
    lines = result.output.splitlines()
    assert any(
        "generators" in line
        and "eos9p57" in line
        and "ok" in line
        and "CReM Growth" in line
        for line in lines
    )
    assert any("similarity" in line and "eos8nyg" in line for line in lines)
    assert any(
        "annotators" in line and "eos4zfy" in line and "MAIP" in line for line in lines
    )


def test_an_unfinished_environment_is_visible(store_root):
    _fetched(store_root, "eos9taz", complete=False)

    result = _catalog()

    assert any(
        "eos9taz" in line and "unfinished" in line
        for line in result.output.splitlines()
    )


def test_an_empty_store_says_how_to_fill_it(store_root):
    result = _catalog()

    assert result.exit_code == 0, result.output
    assert "No models here yet" in result.output
    assert "chemsampler fetch" in result.output
    assert not store_root.exists()  # looking never creates the store


def test_a_project_marks_the_models_it_needs(tmp_path, store_root, monkeypatch):
    create_project(tmp_path / "proj")
    # the project's template wants eos9taz, eos8nyg and eos4zfy; eos9p57 is spare
    _fetched(store_root, "eos9taz", "eos9p57", "eos8nyg", "eos4zfy")
    monkeypatch.chdir(tmp_path / "proj")

    result = _catalog()

    rows = [
        line for line in result.output.splitlines() if re.match(r"^\* \w+ +eos", line)
    ]
    assert {line.split()[2] for line in rows} == {"eos9taz", "eos8nyg", "eos4zfy"}
    assert "* needed by this project" in result.output  # the legend, not a model


def test_models_a_project_needs_but_has_not_fetched_are_called_out(
    tmp_path, store_root, monkeypatch
):
    create_project(tmp_path / "proj")
    _fetched(store_root, "eos9taz")  # eos8nyg and eos4zfy are missing
    monkeypatch.chdir(tmp_path / "proj")

    result = _catalog()

    assert any(
        line.startswith("*") and "eos8nyg" in line
        for line in result.output.splitlines()
    )
    assert "2 model(s) this project needs are not ready" in result.output
    assert "eos8nyg" in result.output and "chemsampler fetch" in result.output


def test_a_cloned_model_whose_environment_failed_counts_as_not_ready(
    tmp_path, store_root, monkeypatch
):
    create_project(tmp_path / "proj")
    _fetched(store_root, "eos9taz", "eos8nyg")
    _fetched(store_root, "eos4zfy", complete=False)
    monkeypatch.chdir(tmp_path / "proj")

    result = _catalog()

    assert "1 model(s) this project needs are not ready (eos4zfy)" in result.output


def test_unreadable_project_inputs_are_explained_not_fatal(
    tmp_path, store_root, monkeypatch
):
    root = tmp_path / "proj"
    create_project(root)
    (root / "input" / "annotators.csv").unlink()
    _fetched(store_root, "eos9taz")
    monkeypatch.chdir(root)

    result = _catalog()

    assert result.exit_code == 0, result.output
    assert "this project's inputs could not be read" in result.output
    assert any("eos9taz" in line for line in result.output.splitlines())


def test_outside_a_project_there_is_nothing_to_warn_about(
    tmp_path, store_root, monkeypatch
):
    _fetched(store_root, "eos9taz")
    monkeypatch.chdir(tmp_path)

    result = _catalog()

    assert "could not be read" not in result.output and "needed" not in result.output


def test_the_role_filter_applies_to_the_local_listing(store_root):
    _fetched(store_root, "eos9taz", "eos8nyg", "eos4zfy")

    result = _catalog("--role", "similarity")

    assert "eos8nyg" in result.output
    assert "eos9taz" not in result.output and "eos4zfy" not in result.output


def test_a_project_is_listed_from_the_folders_it_names(tmp_path):
    clones = tmp_path / "shared" / "models"
    (clones / "eos9taz").mkdir(parents=True)
    env = tmp_path / "shared" / "envs" / "eos9taz"
    env.mkdir(parents=True)
    write_marker(env, model="eos9taz")
    create_project(tmp_path / "proj", models=clones, envs=tmp_path / "shared" / "envs")

    result = _catalog("-p", str(tmp_path / "proj"))

    assert f"Models in {clones}" in result.output
    assert any(
        "eos9taz" in line and "ok" in line for line in result.output.splitlines()
    )


def test_a_folder_named_alone_has_no_environments_to_report(tmp_path):
    clones = tmp_path / "bundles"
    (clones / "eos9taz").mkdir(parents=True)
    create_project(tmp_path / "proj", models=clones)

    result = _catalog("-p", str(tmp_path / "proj"))

    assert "environments in" not in result.output
    assert any("eos9taz" in line for line in result.output.splitlines())


def test_a_broken_catalog_does_not_stop_the_local_listing(store_root, hub_catalog):
    _fetched(store_root, "eos9taz")
    hub_catalog.write_text("garbage")

    result = _catalog()

    assert result.exit_code == 0, result.output
    assert "chemsampler catalog --update" in result.output  # how to repair it
    assert any(
        "eos9taz" in line and "unknown" in line for line in result.output.splitlines()
    )


# --- what the Hub offers ---------------------------------------------------------------


def test_the_hub_view_counts_each_kind_and_what_is_on_this_machine(store_root):
    _fetched(store_root, "eos9taz", "eos8nyg")

    result = _catalog("--hub")

    assert result.exit_code == 0, result.output
    assert "catalog built" in result.output and "ersilia 0.1.62" in result.output
    assert any(
        line.split()[:3] == ["generators", "2", "models,"]
        and "1 on this machine" in line
        for line in result.output.splitlines()
    )
    assert any(
        "similarity" in line and "2 models" in line
        for line in result.output.splitlines()
    )
    assert "--role" in result.output  # how to see the list


def test_the_hub_view_lists_one_kind_and_marks_it(store_root):
    _fetched(store_root, "eos8nyg")

    result = _catalog("--hub", "--role", "similarity")

    lines = result.output.splitlines()
    assert any("eos8nyg" in line and "on this machine" in line for line in lines)
    assert any(
        "eos2a9n" in line and "not supported by chemsampler" in line for line in lines
    )
    assert not any("eos8nyg" in line and "not supported" in line for line in lines)
    assert "2 similarity." in result.output


def test_the_hub_view_needs_a_readable_catalog(hub_catalog):
    hub_catalog.write_text("garbage")

    result = _catalog("--hub")

    assert result.exit_code == 1


# --- refreshing the catalog ------------------------------------------------------------


@pytest.mark.parametrize("flags", [["--hub"], ["--role", "generators"]])
def test_update_does_not_combine_with_listing_options(flags):
    result = _catalog("--update", *flags)

    assert result.exit_code == 2
    assert "cannot be combined" in result.output


def test_update_never_reads_the_store(store_root, monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("refreshing the catalog must not look at the store")

    monkeypatch.setattr(config, "locate_models", fail)
    monkeypatch.setattr(
        catalog,
        "update_catalog",
        lambda: catalog.CatalogUpdate(
            catalog.read_catalog(),
            catalog.read_catalog(),
            catalog.CatalogDiff((), ()),
            False,
        ),
    )

    result = _catalog("--update")

    assert result.exit_code == 0, result.output
    assert "Catalog is up to date" in result.output


def test_a_marker_file_that_names_its_own_models_is_listed_from_there(
    tmp_path, store_root, monkeypatch
):
    # a project made before the store existed keeps its models beside it
    root = tmp_path / "old"
    create_project(root)
    (root / "chemsampler.json").write_text(
        json.dumps({"path_to_model_repositories": "models", "path_to_envs": "envs"})
    )
    (root / "models" / "eos9taz").mkdir(parents=True)
    _fetched(store_root, "eos4zfy")  # in the store, but not this project's
    monkeypatch.chdir(root)

    result = _catalog()

    lines = result.output.splitlines()
    assert f"Models in {root.resolve() / 'models'}" in result.output
    assert any("eos9taz" in line for line in lines)
    # eos4zfy is complete in the store, but this project keeps its own models:
    # the store's copy must not be counted
    assert any("eos4zfy" in line and "missing" in line for line in lines)
    assert not any("eos4zfy" in line and " ok " in line for line in lines)


def test_a_listing_prints_only_the_listing(tmp_path, store_root, monkeypatch, caplog):
    # the real model check, which is what logs "Checking models against the catalog"
    monkeypatch.setattr(config, "validate_models", catalog.validate_models)
    create_project(tmp_path / "proj")
    _fetched(store_root, "eos9taz", "eos8nyg", "eos4zfy")
    monkeypatch.chdir(tmp_path / "proj")

    with caplog.at_level("INFO"):
        result = _catalog()

    assert result.exit_code == 0, result.output
    assert any("eos9taz" in line for line in result.output.splitlines())  # it did look
    assert not [r for r in caplog.records if r.levelname == "INFO"]
    assert "Checking models against the catalog" not in result.output
