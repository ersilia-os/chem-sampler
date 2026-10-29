import json
import re
import subprocess

import pytest

from chemsampler.hub import catalog
from chemsampler.hub.catalog import (
    CatalogEntry,
    build_catalog,
    diff_catalogs,
    read_catalog,
    update_catalog,
    validate_models,
    write_catalog,
)

GENERATOR = CatalogEntry("generators", "eos9taz", "MoLeR", "Sampling", "Generation")
SIMILARITY = CatalogEntry(
    "similarity", "eos8nyg", "Reference library", "Sampling", "Similarity search"
)
OTHER_SIMILARITY = CatalogEntry(
    "similarity", "eos2a9n", "ChEMBL search", "Sampling", "Similarity search"
)
ANNOTATOR = CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", "")
ENTRIES = [GENERATOR, SIMILARITY, OTHER_SIMILARITY, ANNOTATOR]


@pytest.fixture
def catalog_file(tmp_path, monkeypatch):
    """A catalog file in a temporary folder, standing in for the shipped one."""
    path = tmp_path / "catalog.csv"
    write_catalog(ENTRIES, path)
    monkeypatch.setattr(catalog, "CATALOG_PATH", path)
    return path


# --- the file -----------------------------------------------------------------------


def test_write_then_read_round_trips_and_sorts(tmp_path):
    path = tmp_path / "catalog.csv"
    tricky = CatalogEntry(
        "annotators", "eos1aaa", 'A "quoted", title', "Annotation", ""
    )

    write_catalog([ANNOTATOR, tricky, SIMILARITY, GENERATOR], path)
    read = read_catalog(path)

    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", read.built)
    assert [e.identifier for e in read.entries] == [
        "eos9taz",
        "eos8nyg",
        "eos1aaa",
        "eos4zfy",
    ]
    assert tricky in read.entries
    assert not list(tmp_path.glob("*.tmp"))


def test_a_missing_catalog_says_how_to_build_it(tmp_path):
    with pytest.raises(ValueError, match="chemsampler catalog --update"):
        read_catalog(tmp_path / "nope.csv")


@pytest.mark.parametrize(
    "text",
    [
        "role,identifier,title,task,subtask\n",  # no header line
        "# chemsampler catalog, built 2026-10-29, ersilia 0.1.62\nwrong,columns\n",
        (
            "# chemsampler catalog, built 2026-10-29, ersilia 0.1.62\n"
            "role,identifier,title,task,subtask\nboss,eos1,T,Sampling,Generation\n"
        ),
    ],
)
def test_an_invalid_catalog_says_how_to_rebuild_it(tmp_path, text):
    path = tmp_path / "catalog.csv"
    path.write_text(text)

    with pytest.raises(ValueError, match="chemsampler catalog --update"):
        read_catalog(path)


# --- validating models against it ---------------------------------------------------


def test_accepts_models_with_the_right_labels(catalog_file):
    validate_models(
        generators=["eos9taz"], similarity=["eos8nyg"], annotators=["eos4zfy"]
    )


def test_an_unknown_id_names_the_catalog_date_and_the_fix(catalog_file):
    with pytest.raises(ValueError) as error:
        validate_models(generators=["eosXXXX"])

    message = str(error.value)
    assert "eosXXXX (generators): not in the catalog built on" in message
    assert read_catalog(catalog_file).built in message
    assert "chemsampler catalog --update" in message


def test_rejects_models_in_the_wrong_role(catalog_file):
    with pytest.raises(ValueError) as error:
        validate_models(
            generators=["eos8nyg", "eos4zfy"],
            similarity=["eos9taz"],
            annotators=["eos9taz"],
        )

    message = str(error.value)
    assert "eos8nyg (generators)" in message and "put it in similarity" in message
    assert "eos4zfy (generators)" in message and "put it in annotators" in message
    assert "eos9taz (similarity)" in message and "put it in generators" in message
    assert (
        "eos9taz (annotators)" in message
        and "cannot be used as an annotator" in message
    )


def test_rejects_similarity_models_chemsampler_does_not_support(catalog_file):
    with pytest.raises(ValueError) as error:
        validate_models(similarity=["eos2a9n"])

    assert "not meant to be used in chemsampler" in str(error.value)
    assert "Choose from: eos8nyg" in str(error.value)


def test_an_unsupported_similarity_model_in_generators_is_explained(catalog_file):
    with pytest.raises(ValueError, match="similarity model either.*eos8nyg"):
        validate_models(generators=["eos2a9n"])


def test_reports_every_problem_in_one_error(catalog_file):
    with pytest.raises(ValueError) as error:
        validate_models(
            generators=["eosXXXX"], similarity=["eos2a9n"], annotators=["eos9taz"]
        )

    message = str(error.value)
    assert all(model in message for model in ("eosXXXX", "eos2a9n", "eos9taz"))


def test_nothing_to_check_does_not_need_a_catalog(tmp_path, monkeypatch):
    monkeypatch.setattr(catalog, "CATALOG_PATH", tmp_path / "missing.csv")

    validate_models()


# --- building it from the Hub -------------------------------------------------------

SAMPLING_ROWS = [
    {"Identifier": "eos9taz", "Title": "MoLeR", "Task": "Sampling"},
    {"Identifier": "eos8nyg", "Title": "Reference library", "Task": "Sampling"},
    {"Identifier": "eosnew1", "Title": "Future thing", "Task": "Sampling"},
]
CARDS = {
    "eos9taz": {"Subtask": "Generation"},
    "eos8nyg": {"Subtask": "Similarity search"},
    "eosnew1": {"Subtask": "Something else"},
}
ANNOTATION_ROWS = [{"Identifier": "eos4zfy", "Title": "MAIP", "Task": "Annotation"}]


@pytest.fixture
def hub(monkeypatch):
    """A fake Hub behind the single `_run_ersilia` seam."""

    def fake(*args):
        if args[0] == "--more":
            return {"sampling": SAMPLING_ROWS, "annotation": ANNOTATION_ROWS}[args[2]]
        if args[1] not in CARDS:
            raise catalog._NoOutput(
                f"Error: No metadata found for model ID '{args[1]}'"
            )
        return CARDS[args[1]]

    monkeypatch.setattr(catalog, "_run_ersilia", fake)
    return fake


def test_build_gives_each_model_the_role_its_labels_imply(hub):
    built = {e.identifier: e for e in build_catalog()}

    assert (built["eos9taz"].role, built["eos9taz"].subtask) == (
        "generators",
        "Generation",
    )
    assert built["eos8nyg"].role == "similarity"
    assert (built["eos4zfy"].role, built["eos4zfy"].subtask) == ("annotators", "")
    assert "eosnew1" not in built  # a Subtask no role uses


def test_build_fails_when_the_hub_lists_a_model_it_has_no_card_for(hub, monkeypatch):
    monkeypatch.delitem(CARDS, "eos9taz")

    with pytest.raises(RuntimeError, match="has no card"):
        build_catalog()


# --- comparing and updating ---------------------------------------------------------


def test_diff_reports_added_and_removed_models(tmp_path):
    path = tmp_path / "catalog.csv"
    write_catalog([GENERATOR, ANNOTATOR], path)
    old = read_catalog(path)

    diff = diff_catalogs(old, [GENERATOR, SIMILARITY])

    assert diff.added == (SIMILARITY,) and diff.removed == (ANNOTATOR,)
    assert not diff.empty
    assert diff_catalogs(old, [GENERATOR, ANNOTATOR]).empty
    assert diff_catalogs(None, [GENERATOR]).added == (GENERATOR,)


def test_update_leaves_the_file_alone_when_nothing_changed(hub, tmp_path):
    path = tmp_path / "catalog.csv"
    write_catalog(build_catalog(), path)
    before = path.read_bytes()

    result = update_catalog(path)

    assert result.diff.empty and not result.written
    assert path.read_bytes() == before


def test_update_rewrites_the_file_and_reports_what_changed(hub, tmp_path):
    path = tmp_path / "catalog.csv"
    write_catalog([GENERATOR, ANNOTATOR, OTHER_SIMILARITY], path)

    result = update_catalog(path)

    assert result.written
    assert {e.identifier for e in result.diff.added} == {"eos8nyg"}
    assert {e.identifier for e in result.diff.removed} == {"eos2a9n"}
    assert {e.identifier for e in read_catalog(path).entries} == {
        "eos9taz",
        "eos8nyg",
        "eos4zfy",
    }


def test_update_creates_a_missing_file(hub, tmp_path):
    path = tmp_path / "catalog.csv"

    result = update_catalog(path)

    assert result.written and result.previous is None
    assert path.exists()


def test_update_rebuilds_a_corrupt_file(hub, tmp_path):
    path = tmp_path / "catalog.csv"
    path.write_text("garbage")

    result = update_catalog(path)

    assert result.written and result.previous is None
    assert read_catalog(path).count("generators") == 1


def test_update_keeps_the_existing_file_when_the_hub_is_unreachable(
    tmp_path, monkeypatch
):
    path = tmp_path / "catalog.csv"
    write_catalog(ENTRIES, path)
    before = path.read_bytes()

    def unreachable(*args):
        raise catalog._unreachable("connection refused")

    monkeypatch.setattr(catalog, "_run_ersilia", unreachable)

    with pytest.raises(RuntimeError, match="could not query the Ersilia Model Hub"):
        update_catalog(path)

    assert path.read_bytes() == before


# --- the ersilia wrapper itself -----------------------------------------------------


def _fake_subprocess(returncode=0, stdout="", stderr="", payload=None):
    def run(command, **kwargs):
        if payload is not None:
            with open(command[-1], "w") as f:
                json.dump(payload, f)
        return subprocess.CompletedProcess(command, returncode, stdout, stderr)

    return run


@pytest.fixture
def ersilia(monkeypatch):
    monkeypatch.setattr(catalog, "_ersilia_executable", lambda: "ersilia")

    def use(**kwargs):
        monkeypatch.setattr(catalog.subprocess, "run", _fake_subprocess(**kwargs))

    return use


def test_reads_the_list_written_by_ersilia(ersilia):
    ersilia(payload=[{"Identifier": "eos9taz"}, {"Identifier": "eos8nyg"}])

    rows = catalog._list_models("sampling")

    assert [r["Identifier"] for r in rows] == ["eos9taz", "eos8nyg"]


def test_unknown_id_exits_zero_without_a_file(ersilia):
    ersilia(stderr="Error: No metadata found for model ID 'eosXXXX'")

    assert catalog._card("eosXXXX") is None


def test_missing_file_without_the_unknown_id_text_means_unreachable(ersilia):
    ersilia(stdout="Error fetching model metadata: connection refused")

    with pytest.raises(RuntimeError, match="could not query the Ersilia Model Hub"):
        catalog._card("eos9taz")


def test_a_list_that_was_not_written_is_an_error(ersilia):
    ersilia()

    with pytest.raises(RuntimeError, match="could not query the Ersilia Model Hub"):
        catalog._list_models("sampling")


def test_a_failing_ersilia_is_an_error(ersilia):
    ersilia(returncode=1, stderr="boom")

    with pytest.raises(RuntimeError, match="exited with 1: boom"):
        catalog._list_models("sampling")


def test_a_timeout_is_an_error(monkeypatch):
    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, 60)

    monkeypatch.setattr(catalog, "_ersilia_executable", lambda: "ersilia")
    monkeypatch.setattr(catalog.subprocess, "run", timeout)

    with pytest.raises(RuntimeError, match="timed out"):
        catalog._list_models("annotation")
