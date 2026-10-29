import json

import pytest
from click.testing import CliRunner

from chemsampler.cli import cli
from chemsampler.config import load_inputs
from chemsampler.workspace import (
    TEMPLATES,
    WORKSPACE_FILE,
    Workspace,
    create_workspace,
    find_workspace,
)


def test_setup_creates_the_layout(tmp_path):
    created = create_workspace(tmp_path / "run")

    root = tmp_path / "run"
    for folder in ("input", "models", "envs", "rounds", "logs"):
        assert (root / folder).is_dir()
    assert (root / WORKSPACE_FILE).is_file()
    assert (root / ".gitignore").read_text().split() == [
        "models/",
        "envs/",
        "rounds/",
        "logs/",
    ]
    for name in TEMPLATES:
        assert (root / "input" / name).is_file()
    assert WORKSPACE_FILE in created and "input/" in created


def test_setup_is_safe_to_repeat_and_keeps_edited_inputs(tmp_path):
    create_workspace(tmp_path)
    (tmp_path / "input" / "generators.csv").write_text("generators\neos8vud\n")

    again = create_workspace(tmp_path)

    assert again == []
    assert (
        tmp_path / "input" / "generators.csv"
    ).read_text() == "generators\neos8vud\n"


def test_setup_fills_in_only_what_is_missing(tmp_path):
    create_workspace(tmp_path)
    (tmp_path / "envs").rmdir()
    (tmp_path / "input" / "similarity.csv").unlink()

    assert sorted(create_workspace(tmp_path)) == ["envs/", "input/similarity.csv"]


def test_a_custom_layout_in_the_marker_file_is_respected(tmp_path):
    (tmp_path / WORKSPACE_FILE).write_text(json.dumps({"path_to_models": "csvs"}))

    create_workspace(tmp_path)
    workspace = Workspace.load(tmp_path)

    assert workspace.input == tmp_path / "csvs"
    assert (tmp_path / "csvs" / "annotators.csv").is_file()
    assert workspace.models == tmp_path / "models"  # the default for a key it omits


def test_the_workspace_is_found_from_any_subfolder(tmp_path):
    create_workspace(tmp_path)
    deep = tmp_path / "models" / "eos9taz" / "model"
    deep.mkdir(parents=True)

    assert find_workspace(deep) == tmp_path.resolve()
    assert find_workspace(tmp_path / "input") == tmp_path.resolve()


def test_no_workspace_is_found_outside_one(tmp_path):
    assert find_workspace(tmp_path) is None


def test_a_broken_marker_file_is_explained(tmp_path):
    (tmp_path / WORKSPACE_FILE).write_text("[1, 2]")
    with pytest.raises(ValueError, match="expected a JSON object"):
        Workspace.load(tmp_path)

    (tmp_path / WORKSPACE_FILE).write_text("{oops")
    with pytest.raises(ValueError, match="not valid JSON"):
        Workspace.load(tmp_path)


def test_the_templates_are_valid_inputs_for_the_shipped_catalog(tmp_path):
    create_workspace(tmp_path)

    inputs = load_inputs(path_to_models=str(tmp_path / "input"), backend="ersilia")

    assert inputs.similarity.model_id == "eos8nyg"
    assert [g.model_id for g in inputs.generators.generators] == ["eos9taz"]
    assert [a.annotator_id for a in inputs.annotators] == ["eos4zfy"]


def test_the_setup_command_reports_what_it_made_and_what_it_found(tmp_path):
    target = tmp_path / "run"

    first = CliRunner().invoke(cli, ["setup", str(target)])
    second = CliRunner().invoke(cli, ["setup", str(target)])

    assert first.exit_code == 0, first.output
    assert "Workspace ready" in first.output and "input/" in first.output
    assert "not move or rename" in first.output
    assert second.exit_code == 0, second.output
    assert "already a workspace" in second.output


def test_the_setup_command_works_in_the_current_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(cli, ["setup"])

    assert result.exit_code == 0, result.output
    assert (tmp_path / WORKSPACE_FILE).is_file()
