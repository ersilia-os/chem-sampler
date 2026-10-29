import json

import pytest
from click.testing import CliRunner

from chemsampler.cli import cli
from chemsampler.config import load_inputs
from chemsampler.project import (
    PROJECT_FILE,
    TEMPLATES,
    Project,
    create_project,
    find_project,
)


def test_init_creates_only_what_a_project_owns(tmp_path):
    created = create_project(tmp_path / "run")

    root = tmp_path / "run"
    assert sorted(p.name for p in root.iterdir()) == [
        ".gitignore",
        PROJECT_FILE,
        "input",
        "rounds",
    ]
    assert json.loads((root / PROJECT_FILE).read_text()) == {
        "path_to_models": "input",
        "path_to_rounds": "rounds",
    }
    assert (root / ".gitignore").read_text().split() == ["rounds/"]
    for name in TEMPLATES:
        assert (root / "input" / name).is_file()
    assert PROJECT_FILE in created and "input/" in created and "rounds/" in created


def test_a_project_has_no_models_or_envs_of_its_own(tmp_path):
    # those are the store's: expensive to build, so shared between projects
    create_project(tmp_path)

    project = Project.load(tmp_path)

    assert project.models is None and project.envs is None
    assert not (tmp_path / "models").exists() and not (tmp_path / "envs").exists()


def test_init_is_safe_to_repeat_and_keeps_edited_inputs(tmp_path):
    create_project(tmp_path)
    (tmp_path / "input" / "generators.csv").write_text("generators\neos8vud\n")

    again = create_project(tmp_path)

    assert again == []
    assert (
        tmp_path / "input" / "generators.csv"
    ).read_text() == "generators\neos8vud\n"


def test_init_fills_in_only_what_is_missing(tmp_path):
    create_project(tmp_path)
    (tmp_path / "rounds").rmdir()
    (tmp_path / "input" / "similarity.csv").unlink()

    assert sorted(create_project(tmp_path)) == ["input/similarity.csv", "rounds/"]


def test_a_custom_layout_in_the_marker_file_is_respected(tmp_path):
    (tmp_path / PROJECT_FILE).write_text(json.dumps({"path_to_models": "csvs"}))

    create_project(tmp_path)
    project = Project.load(tmp_path)

    assert project.input == tmp_path / "csvs"
    assert (tmp_path / "csvs" / "annotators.csv").is_file()
    assert project.rounds == tmp_path / "rounds"  # the default for a key it omits


def test_a_marker_that_names_its_own_models_keeps_them_apart_from_the_store(tmp_path):
    # how a project made before the store existed keeps working
    (tmp_path / PROJECT_FILE).write_text(
        json.dumps({"path_to_model_repositories": "models", "path_to_envs": "envs"})
    )

    project = Project.load(tmp_path)

    assert project.models == tmp_path / "models"
    assert project.envs == tmp_path / "envs"


def test_the_project_is_found_from_any_subfolder(tmp_path):
    create_project(tmp_path)
    deep = tmp_path / "rounds" / "2026-10-29_16-51-01"
    deep.mkdir(parents=True)

    assert find_project(deep) == tmp_path.resolve()
    assert find_project(tmp_path / "input") == tmp_path.resolve()


def test_no_project_is_found_outside_one(tmp_path):
    assert find_project(tmp_path) is None


def test_a_broken_marker_file_is_explained(tmp_path):
    (tmp_path / PROJECT_FILE).write_text("[1, 2]")
    with pytest.raises(ValueError, match="expected a JSON object"):
        Project.load(tmp_path)

    (tmp_path / PROJECT_FILE).write_text("{oops")
    with pytest.raises(ValueError, match="not valid JSON"):
        Project.load(tmp_path)


def test_the_templates_are_valid_inputs_for_the_shipped_catalog(tmp_path):
    create_project(tmp_path)

    inputs = load_inputs(path_to_models=str(tmp_path / "input"), backend="ersilia")

    assert inputs.similarity.model_id == "eos8nyg"
    assert [g.model_id for g in inputs.generators.generators] == ["eos9taz"]
    assert [a.annotator_id for a in inputs.annotators] == ["eos4zfy"]


def test_the_init_command_reports_what_it_made_and_what_it_found(tmp_path):
    target = tmp_path / "run"

    first = CliRunner().invoke(cli, ["init", str(target)])
    second = CliRunner().invoke(cli, ["init", str(target)])

    assert first.exit_code == 0, first.output
    assert "Project ready" in first.output and "input/" in first.output
    assert "store" in first.output  # says where the models go
    assert second.exit_code == 0, second.output
    assert "already a project" in second.output


def test_the_init_command_works_in_the_current_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(cli, ["init"])

    assert result.exit_code == 0, result.output
    assert (tmp_path / PROJECT_FILE).is_file()


def test_setup_is_gone(tmp_path):
    result = CliRunner().invoke(cli, ["setup", str(tmp_path / "run")])

    assert result.exit_code == 2 and "No such command" in result.output
