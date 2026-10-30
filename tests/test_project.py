import json
from pathlib import Path

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


# --- the init command ---------------------------------------------------------------


@pytest.fixture
def here(tmp_path, monkeypatch):
    """An empty working directory, apart from the store the suite points at."""
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    return work


def init(*args):
    return CliRunner().invoke(cli, ["init", *args])


def snapshot(folder):
    return {
        p.relative_to(folder): p.read_bytes() if p.is_file() else None
        for p in folder.rglob("*")
    }


def test_init_without_a_path_makes_a_dummy_project_and_nothing_else(here, store_root):
    result = init()

    assert result.exit_code == 0, result.output
    assert [p.name for p in here.iterdir()] == ["chemsampler-test"]
    assert (here / "chemsampler-test" / PROJECT_FILE).is_file()
    assert "Project ready" in result.output and "input/" in result.output
    assert "dummy project" in result.output
    assert "Next: cd chemsampler-test," in result.output
    assert not store_root.exists()  # only fetch makes the store


@pytest.mark.parametrize("flag", ["-p", "--path"])
def test_init_makes_the_folder_it_is_given(here, flag):
    result = init(flag, "a/b/run")

    assert result.exit_code == 0, result.output
    assert (here / "a" / "b" / "run" / PROJECT_FILE).is_file()
    assert "dummy" not in result.output
    assert "Next: cd a/b/run," in result.output


@pytest.mark.parametrize(
    "args, name",
    [([], "chemsampler-test"), (["-p", "mine"], "mine"), (["--path", "."], ".")],
)
def test_init_never_writes_into_a_project_that_exists(here, args, name):
    create_project(here / name)
    (here / name / "input" / "annotators.csv").write_text("edited\n")
    (here / name / "input" / "generators.csv").unlink()
    before = snapshot(here)

    result = init(*args)

    assert result.exit_code == 1
    assert "already exists" in result.stderr
    assert snapshot(here) == before  # not even the missing template comes back


@pytest.mark.parametrize(
    "make",
    [Path.mkdir, Path.touch, lambda path: path.symlink_to(path.parent / "nowhere")],
    ids=["empty folder", "file", "dangling link"],
)
def test_any_path_that_exists_counts_as_taken(here, make):
    make(here / "chemsampler-test")
    before = snapshot(here)

    default = init()
    explicit = init("-p", "chemsampler-test")

    for result in (default, explicit):
        assert result.exit_code == 1
        assert "already exists" in result.stderr
    assert snapshot(here) == before
    assert not (here / "nowhere").exists()  # a link is never followed


def test_the_old_positional_form_is_refused(here):
    result = init("run")

    assert result.exit_code == 2
    assert list(here.iterdir()) == []


def test_a_failed_init_leaves_nothing_that_blocks_a_retry(here, monkeypatch):
    from chemsampler import project

    real = project.create_project
    calls = []

    def fails_once(root):
        calls.append(root)
        if len(calls) == 1:
            (Path(root) / PROJECT_FILE).write_text("{}")
            raise OSError("disk full")
        return real(root)

    monkeypatch.setattr(project, "create_project", fails_once)
    (here / "kept").mkdir()

    failed = init("-p", "kept/run")
    retried = init("-p", "kept/run")

    assert failed.exit_code == 1 and "disk full" in failed.stderr
    assert retried.exit_code == 0, retried.output
    assert len(calls) == 2  # the retry was not refused as "already exists"
    assert (here / "kept" / "run" / "input").is_dir()


def test_a_failed_init_removes_the_half_made_folder(here, monkeypatch):
    from chemsampler import project

    def half_done(root):
        (Path(root) / PROJECT_FILE).write_text("{}")
        raise OSError("disk full")

    monkeypatch.setattr(project, "create_project", half_done)
    (here / "kept").mkdir()

    result = init("-p", "kept/run")

    assert result.exit_code == 1
    assert [p.name for p in here.iterdir()] == ["kept"]
    assert list((here / "kept").iterdir()) == []  # a parent that was there stays


def test_the_next_step_names_the_folder_to_cd_into(here, tmp_path):
    outside = init("-p", str(tmp_path / "elsewhere"))
    spaced = init("-p", "my run")

    assert f"Next: cd {tmp_path / 'elsewhere'}," in outside.output
    assert "Next: cd 'my run'," in spaced.output


def test_init_says_where_the_models_go(here, store_root):
    before = init("-p", "one")
    store_root.mkdir()
    after = init("-p", "two")

    line = f"Models go to the store shared by every project: {store_root}"
    assert f"{line} (`chemsampler fetch` creates it)\n" in before.output
    assert f"{line}\n" in after.output


def test_a_broken_store_setting_does_not_stop_init(here, store_root):
    store_root.mkdir()
    (store_root / "config.json").write_text("[1, 2]")

    result = init("-p", "run")

    assert result.exit_code == 0, result.output
    assert "store settings are not valid" in result.output
    assert (here / "run" / PROJECT_FILE).is_file()


def test_init_notes_a_project_inside_another(here):
    create_project(here / "outer")

    inside = init("-p", "outer/rounds/inner")
    alone = init("-p", "alone")

    assert inside.exit_code == 0, inside.output
    assert f"inside the project {(here / 'outer').resolve()}" in inside.output
    assert "inside the project" not in alone.output


def test_setup_is_gone(tmp_path):
    result = CliRunner().invoke(cli, ["setup", str(tmp_path / "run")])

    assert result.exit_code == 2 and "No such command" in result.output
