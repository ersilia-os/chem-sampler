import json
import subprocess
from pathlib import Path

import pytest

from chemsampler.hub import fetch
from chemsampler.hub.envs import env_path, is_complete, write_marker
from chemsampler.hub.fetch import (
    FetchPlan,
    check_env,
    check_prerequisites,
    fetch_models,
    read_recipe,
    summarize,
)

YML = (
    'python: "3.12"\n'
    "commands:\n"
    '  - ["pip", "torch", "2.9.1+cpu", "--index-url", "https://download.pytorch.org/whl/cpu"]\n'
    '  - ["pip", "numpy", "2.5.3"]\n'
)
DOCKERFILE = (
    "FROM bentoml/model-server:0.11.0-py312\n"
    "RUN pip install torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu\n"
    "RUN pip install transformers==4.57.6 \\\n    huggingface-hub==0.36.2\n"
    "RUN python -m pip install --upgrade pip\n"
)


# --- a fake git, eosvc, conda and pip ---------------------------------------------------


def _done(stdout="", stderr="", code=0):
    return subprocess.CompletedProcess([], code, stdout, stderr)


class FakeSystem:
    """Stands in for every command `fetch` runs, doing their file-system effects."""

    def __init__(self, recipes, access=()):
        self.recipes = recipes  # model id -> (file name, text)
        self.access = set(access)  # ids whose clone carries an access.json
        self.behind, self.ahead, self.dirty = {}, {}, set()
        self.fail = {}  # a text in the command line -> the stderr it fails with
        self.installed, self.conda_installed, self.python = {}, {}, {}
        self.calls = []

    def commands(self, *start):
        return [c for c, _, _ in self.calls if c[: len(start)] == list(start)]

    def __call__(self, argv, *, cwd=None, env=None):
        argv = [str(a) for a in argv]
        self.calls.append((argv, cwd, env))
        line = " ".join(argv)
        for needle, message in self.fail.items():
            if needle in line:
                return _done(stderr=message, code=1)
        if argv[0] == "git":
            return self._git(argv)
        if argv[0] == "conda":
            return self._conda(argv)
        if argv[1:3] == ["-m", "pip"]:
            return self._pip(Path(argv[0]).parent.parent, argv[3:])
        if argv[1] == "-c":
            return _done(self.python.get(Path(argv[0]).parent.parent, "3.12") + "\n")
        return _done()  # eosvc

    def _git(self, argv):
        if argv[1] == "clone":
            dest = Path(argv[-1])
            (dest / ".git").mkdir(parents=True)
            name, text = self.recipes[dest.name]
            (dest / name).write_text(text)
            if dest.name in self.access:
                (dest / "access.json").write_text("{}")
            return _done()
        model = Path(argv[2]).name
        return {
            "fetch": lambda: _done(),
            "rev-list": lambda: _done(
                f"{self.ahead.get(model, 0)}\t{self.behind.get(model, 0)}\n"
            ),
            "status": lambda: _done(" M code.py\n" if model in self.dirty else ""),
            "rev-parse": lambda: _done("abc1234\n"),
            "pull": lambda: self._pull(model),
        }[argv[3]]()

    def _pull(self, model):
        self.behind[model] = 0
        return _done()

    def _conda(self, argv):
        if argv[1] == "create":
            prefix = Path(argv[3])
            (prefix / "bin").mkdir(parents=True)
            self.installed[prefix], self.conda_installed[prefix] = {}, {}
            self.python[prefix] = argv[4].split("=")[1]
        elif argv[1] == "list":
            installed = self.conda_installed.get(Path(argv[3]), {})
            return _done(
                json.dumps([{"name": n, "version": v} for n, v in installed.items()])
            )
        elif argv[1] == "install":
            prefix = Path(argv[3])
            for token in argv[4:]:
                if "=" in token and not token.startswith("-"):
                    name, version = token.split("=")
                    self.conda_installed[prefix][name] = version
        return _done()

    def _pip(self, prefix, args):
        if args[0] == "list":
            installed = self.installed.get(prefix, {})
            return _done(
                json.dumps([{"name": n, "version": v} for n, v in installed.items()])
            )
        if args[0] == "install":
            for token in args[1:]:
                if "==" in token:
                    name, version = token.split("==")
                    self.installed[prefix][name] = version
        return _done()


@pytest.fixture
def system(monkeypatch):
    fake = FakeSystem(
        {"eos1aaa": ("install.yml", YML), "eos2bbb": ("Dockerfile", DOCKERFILE)}
    )
    monkeypatch.setattr(fetch, "_run", fake)
    monkeypatch.setattr(fetch, "_eosvc", lambda: "eosvc")
    monkeypatch.setattr(fetch, "free_disk_gb", lambda path: 500.0)
    return fake


@pytest.fixture
def plan(tmp_path):
    return FetchPlan(
        model_ids=("eos1aaa",),
        models_dir=tmp_path / "models",
        envs_dir=tmp_path / "envs",
        logs_dir=tmp_path / "logs" / "fetch",
    )


def _say(messages):
    return messages.append


# --- recipes -----------------------------------------------------------------------------


def test_an_install_yml_gives_steps_and_pins(tmp_path):
    (tmp_path / "install.yml").write_text(YML)

    recipe = read_recipe(tmp_path)

    assert recipe.python == "3.12" and recipe.source == "install.yml"
    assert recipe.pip_pins == {"torch": "2.9.1+cpu", "numpy": "2.5.3"}
    assert recipe.steps[0] == (
        "pip",
        [
            "install",
            "torch==2.9.1+cpu",
            "--index-url",
            "https://download.pytorch.org/whl/cpu",
        ],
    )
    assert recipe.unchecked == 0


def test_lines_that_cannot_be_compared_are_counted_not_pinned(tmp_path):
    (tmp_path / "install.yml").write_text(
        'python: "3.10"\ncommands:\n'
        '  - ["pip", "git+https://github.com/x/y.git@abc"]\n'
        '  - ["pip", "unpinned"]\n'
        '  - ["conda", "rdkit", "2023.9.4", "conda-forge", "-y"]\n'
        '  - ["conda", "plain"]\n'
    )

    recipe = read_recipe(tmp_path)

    assert recipe.pip_pins == {} and recipe.conda_pins == {"rdkit": "2023.9.4"}
    assert recipe.unchecked == 3
    assert ("conda", ["-c", "conda-forge", "rdkit=2023.9.4"]) in recipe.steps


def test_a_dockerfile_gives_the_python_version_and_pins(tmp_path):
    (tmp_path / "Dockerfile").write_text(DOCKERFILE)

    recipe = read_recipe(tmp_path)

    assert recipe.python == "3.12" and recipe.source == "Dockerfile"
    assert recipe.pip_pins == {
        "torch": "2.9.1",
        "transformers": "4.57.6",
        "huggingface-hub": "0.36.2",
    }
    assert [kind for kind, _ in recipe.steps] == ["pip", "pip", "pip"]
    assert recipe.unchecked == 1  # `pip install --upgrade pip` has nothing to pin


def test_dockerfile_shell_lines_run_as_shell_and_conda_gets_yes(tmp_path):
    (tmp_path / "Dockerfile").write_text(
        "FROM bentoml/model-server:0.11.0-py310\n"
        "RUN curl -fsSL https://example.org/install.sh | bash\n"
        "RUN conda install -c conda-forge opam\n"
    )

    recipe = read_recipe(tmp_path)

    assert recipe.steps == (
        ("shell", "curl -fsSL https://example.org/install.sh | bash"),
        ("shell", "conda install -c conda-forge opam -y"),
    )
    assert recipe.unchecked == 2


def test_a_model_without_a_recipe_or_with_a_broken_one_is_explained(tmp_path):
    with pytest.raises(ValueError, match="no install.yml or Dockerfile"):
        read_recipe(tmp_path)
    (tmp_path / "install.yml").write_text("python: 3.12\n")
    with pytest.raises(ValueError, match="cannot read the recipe"):
        read_recipe(tmp_path)
    (tmp_path / "install.yml").unlink()
    (tmp_path / "Dockerfile").write_text("FROM ubuntu\n")
    with pytest.raises(ValueError, match="python version"):
        read_recipe(tmp_path)


# --- a first fetch --------------------------------------------------------------------------


def test_a_first_fetch_clones_downloads_and_builds(system, plan):
    system.access = {"eos1aaa"}
    messages = []

    [result] = fetch_models(plan, say=_say(messages))

    model_dir, env = plan.models_dir / "eos1aaa", plan.envs_dir / "eos1aaa"
    assert result.ok and (result.clone, result.checkpoints, result.env) == (
        "cloned",
        "downloaded",
        "built",
    )
    assert system.commands("git", "clone")[0][-2:] == [
        "https://github.com/ersilia-os/eos1aaa.git",
        str(model_dir),
    ]
    [eosvc_call] = [c for c in system.calls if c[0][0] == "eosvc"]
    assert eosvc_call[0] == ["eosvc", "download", "--path", "."]
    assert (
        eosvc_call[1] == model_dir
    )  # eosvc reads the model id from this folder's name
    assert system.commands("conda", "create")[0][2:] == [
        "-p",
        str(env),
        "python=3.12",
        "-y",
    ]
    assert is_complete(env)
    assert (
        json.loads((env / ".chemsampler-env.json").read_text())["commit"] == "abc1234"
    )
    assert "$ git clone" in (plan.logs_dir / "eos1aaa.log").read_text()


def test_pip_always_runs_through_the_environments_own_python(system, plan):
    fetch_models(plan, say=_say([]))

    pip_calls = [
        c for c, _, _ in system.calls if c[1:3] == ["-m", "pip"] and c[3] == "install"
    ]
    assert len(pip_calls) == 2
    assert all(
        c[0] == str(plan.envs_dir / "eos1aaa" / "bin" / "python") for c in pip_calls
    )
    assert not [c for c, _, _ in system.calls if c[0] == "pip"]
    assert all(
        env["PYTHONNOUSERSITE"] == "1"
        for c, _, env in system.calls
        if c[1:3] == ["-m", "pip"]
    )


def test_a_model_without_access_json_needs_no_download(system, plan):
    [result] = fetch_models(plan, say=_say([]))

    assert result.checkpoints == "none needed"
    assert not [c for c, _, _ in system.calls if c[0] == "eosvc"]


def test_a_dockerfile_model_builds_its_pip_steps_too(system, plan):
    plan = FetchPlan(("eos2bbb",), plan.models_dir, plan.envs_dir, plan.logs_dir)

    [result] = fetch_models(plan, say=_say([]))

    assert result.env == "built"
    assert system.installed[plan.envs_dir / "eos2bbb"]["transformers"] == "4.57.6"


# --- running it again -------------------------------------------------------------------------


def test_fetching_again_changes_nothing(system, plan):
    system.access = {"eos1aaa"}
    fetch_models(plan, say=_say([]))
    system.calls.clear()

    [result] = fetch_models(plan, say=_say([]))

    assert (result.clone, result.checkpoints, result.env) == (
        "up to date",
        "present",
        "ok",
    )
    assert not system.commands("git", "clone") and not system.commands("conda")
    assert not [c for c, _, _ in system.calls if c[0] == "eosvc"]
    assert not [c for c in system.commands("git") if "pull" in c]


def test_a_clone_that_is_behind_is_reported_and_left_alone(system, plan):
    fetch_models(plan, say=_say([]))
    system.behind["eos1aaa"] = 3

    [result] = fetch_models(plan, say=_say([]))

    assert result.clone == "3 behind"
    assert not [c for c in system.commands("git") if "pull" in c]


def test_update_pulls_a_clone_that_is_behind_fast_forward_only(system, plan):
    fetch_models(plan, say=_say([]))
    system.behind["eos1aaa"] = 3

    [result] = fetch_models(plan, update=True, say=_say([]))

    assert result.clone == "updated"
    [pull] = [c for c in system.commands("git") if "pull" in c]
    assert "--ff-only" in pull


@pytest.mark.parametrize("trouble", ["dirty", "ahead", "diverged"])
def test_update_never_touches_a_clone_with_local_work(system, plan, trouble):
    fetch_models(plan, say=_say([]))
    system.behind["eos1aaa"] = 2
    if trouble == "dirty":
        system.dirty.add("eos1aaa")
    else:
        system.ahead["eos1aaa"] = 1
    if trouble == "ahead":
        system.behind["eos1aaa"] = 0

    [result] = fetch_models(plan, update=True, say=_say([]))

    assert not [c for c in system.commands("git") if "pull" in c]
    assert result.notes == ["not updated: local changes or local commits"]


def test_a_folder_that_is_not_a_clone_is_not_overwritten(system, plan):
    (plan.models_dir / "eos1aaa").mkdir(parents=True)
    (plan.models_dir / "eos1aaa" / "mine.txt").write_text("keep")

    [result] = fetch_models(plan, say=_say([]))

    assert not result.ok and "not a git clone" in result.error
    assert (plan.models_dir / "eos1aaa" / "mine.txt").read_text() == "keep"


# --- environments that already exist ----------------------------------------------------------


def _built_env(system, plan, installed, python="3.12", model="eos1aaa"):
    env = env_path(plan.envs_dir, model)
    (env / "bin").mkdir(parents=True)
    write_marker(env, model=model)
    system.installed[env], system.python[env] = installed, python
    return env


def _with_clone(system, plan, model="eos1aaa"):
    """Give `model` a clone, without building anything."""
    system(["git", "clone", "x", str(plan.models_dir / model)])


def test_an_environment_that_differs_is_reported_not_rebuilt(system, plan):
    _with_clone(system, plan)
    env = _built_env(system, plan, {"torch": "2.8.0", "numpy": "2.5.3"})
    system.calls.clear()

    [result] = fetch_models(plan, say=_say([]))

    assert result.ok and result.env == "differs"
    assert result.differences == ("torch: installed 2.8.0, recipe 2.9.1+cpu",)
    assert not system.commands("conda") and env.exists()
    report = "\n".join(summarize([result]))
    assert "chemsampler fetch --rebuild" in report


@pytest.mark.parametrize(
    "installed, pinned, matches",
    [
        ("2.9.1", "2.9.1", True),
        (
            "2.9.1+cpu",
            "2.9.1",
            True,
        ),  # the CPU index serves a local label for a plain pin
        ("2.9.1+cu128", "2.9.1", True),
        ("2.9.1+cpu", "2.9.1+cpu", True),
        ("2.9.1", "2.9.1+cpu", False),  # a pin that names a label needs it
        ("2.9.1+cu128", "2.9.1+cpu", False),
        ("2.8.0", "2.9.1", False),
        ("2.9.1.post1", "2.9.1", False),
    ],
)
def test_versions_are_compared_the_way_pip_judges_them(installed, pinned, matches):
    assert fetch._version_matches(installed, pinned) is matches


def test_a_cpu_wheel_satisfies_a_plain_pin(system, plan):
    plan = FetchPlan(("eos2bbb",), plan.models_dir, plan.envs_dir, plan.logs_dir)
    _with_clone(system, plan, "eos2bbb")
    _built_env(
        system,
        plan,
        # the Dockerfile pins torch==2.9.1 and installs it from the CPU index
        {"torch": "2.9.1+cpu", "transformers": "4.57.6", "huggingface-hub": "0.36.2"},
        model="eos2bbb",
    )

    [result] = fetch_models(plan, say=_say([]))

    assert result.env == "ok" and result.differences == ()


def test_a_missing_package_and_a_wrong_python_are_differences(system, plan):
    _with_clone(system, plan)
    _built_env(system, plan, {"torch": "2.9.1+cpu"}, python="3.10")

    [result] = fetch_models(plan, say=_say([]))

    assert result.differences == (
        "python: installed 3.10, recipe 3.12",
        "numpy: not installed, recipe 2.5.3",
    )


def test_package_names_are_compared_in_their_normalized_form(system, plan):
    plan = FetchPlan(("eos2bbb",), plan.models_dir, plan.envs_dir, plan.logs_dir)
    _with_clone(system, plan, "eos2bbb")
    _built_env(
        system,
        plan,
        {"torch": "2.9.1", "transformers": "4.57.6", "huggingface_hub": "0.36.2"},
        model="eos2bbb",
    )

    [result] = fetch_models(plan, say=_say([]))

    assert result.env == "ok"  # huggingface_hub is huggingface-hub


def test_rebuild_rebuilds_only_the_environments_that_differ(system, plan):
    plan = FetchPlan(
        ("eos1aaa", "eos2bbb"), plan.models_dir, plan.envs_dir, plan.logs_dir
    )
    for model in plan.model_ids:
        _with_clone(system, plan, model)
    _built_env(system, plan, {"torch": "2.8.0", "numpy": "2.5.3"})
    _built_env(
        system,
        plan,
        {"torch": "2.9.1", "transformers": "4.57.6", "huggingface-hub": "0.36.2"},
        model="eos2bbb",
    )
    system.calls.clear()

    results = fetch_models(plan, rebuild=True, say=_say([]))

    assert [r.env for r in results] == ["rebuilt", "ok"]
    [create] = system.commands("conda", "create")
    assert str(plan.envs_dir / "eos1aaa") in create
    assert system.installed[plan.envs_dir / "eos1aaa"]["torch"] == "2.9.1+cpu"


def test_conda_pins_are_checked_with_one_conda_list(system, plan, tmp_path):
    _with_clone(system, plan)
    (plan.models_dir / "eos1aaa" / "install.yml").write_text(
        'python: "3.12"\ncommands:\n  - ["conda", "rdkit", "2023.9.4", "conda-forge"]\n'
    )
    env = _built_env(system, plan, {})
    system.conda_installed[env] = {"rdkit": "2024.3.1"}

    [result] = fetch_models(plan, say=_say([]))

    assert result.differences == ("rdkit (conda): installed 2024.3.1, recipe 2023.9.4",)
    assert len(system.commands("conda", "list")) == 1


def test_an_unfinished_environment_is_removed_and_rebuilt(system, plan):
    _with_clone(system, plan)
    env = env_path(plan.envs_dir, "eos1aaa")
    env.mkdir(parents=True)
    (env / "half-installed").write_text("x")
    messages = []

    [result] = fetch_models(plan, say=_say(messages))

    assert result.env == "built" and is_complete(env)
    assert not (env / "half-installed").exists()
    assert any("unfinished environment" in m for m in messages)


def test_check_env_reports_a_broken_python(system, plan):
    recipe_dir = plan.models_dir / "eos1aaa"
    recipe_dir.mkdir(parents=True)
    (recipe_dir / "install.yml").write_text(YML)
    env = _built_env(system, plan, {})
    system.fail["-c"] = "no such file"

    check = check_env(env, read_recipe(recipe_dir), fetch._Log(plan.logs_dir / "x.log"))

    assert not check.matches and "does not run" in check.differences[0]


# --- failures ---------------------------------------------------------------------------------


def test_a_failed_step_leaves_no_marker_and_the_next_fetch_rebuilds(system, plan):
    system.fail["numpy==2.5.3"] = "no matching distribution"

    [result] = fetch_models(plan, say=_say([]))

    env = plan.envs_dir / "eos1aaa"
    assert (
        not result.ok
        and "step 2/2" in result.error
        and "no matching distribution" in result.error
    )
    assert env.exists() and not is_complete(env)

    system.fail.clear()
    [again] = fetch_models(plan, say=_say([]))
    assert again.ok and again.env == "built" and is_complete(env)


def test_one_failing_model_does_not_stop_the_others(system, plan):
    plan = FetchPlan(
        ("eos1aaa", "eos2bbb"), plan.models_dir, plan.envs_dir, plan.logs_dir
    )
    system.fail["eos1aaa.git"] = "repository not found"

    first, second = fetch_models(plan, say=_say([]))

    assert not first.ok and "git clone failed" in first.error
    assert second.ok and second.env == "built"


def test_a_failed_checkpoint_download_is_retried_next_time(system, plan):
    system.access = {"eos1aaa"}
    system.fail["eosvc"] = "S3 unreachable"

    [failed] = fetch_models(plan, say=_say([]))
    assert not failed.ok and "eosvc download failed" in failed.error

    system.fail.clear()
    [again] = fetch_models(plan, say=_say([]))
    assert again.ok and again.checkpoints == "downloaded"


def test_conda_create_failures_show_conda_s_own_message(system, plan):
    system.fail["conda create"] = "CondaToSNonInteractiveError: terms of service"

    [result] = fetch_models(plan, say=_say([]))

    assert "CondaToSNonInteractiveError" in result.error


def test_a_low_disk_warns_first(system, plan, monkeypatch):
    monkeypatch.setattr(fetch, "free_disk_gb", lambda path: 3.2)
    messages = []

    fetch_models(plan, say=_say(messages))

    assert messages[0].startswith("warning: only 3.2 GB free")


# --- the report and the prerequisites ------------------------------------------------------------


def test_the_summary_is_a_table_then_what_needs_attention(system, plan):
    _with_clone(system, plan)
    _built_env(system, plan, {"torch": "2.8.0", "numpy": "2.5.3"})
    failed = fetch.ModelResult("eos9zzz", error="git clone failed: not found")

    lines = summarize([*fetch_models(plan, say=_say([])), failed])

    assert lines[0].split() == ["model", "clone", "checkpoints", "environment"]
    assert any(line.startswith("eos1aaa") and "differs" in line for line in lines)
    assert any(line.startswith("eos9zzz") and "FAILED" in line for line in lines)
    assert (
        "eos1aaa: environment differs - torch: installed 2.8.0, recipe 2.9.1+cpu"
        in lines
    )
    assert "eos9zzz: git clone failed: not found" in lines


def test_missing_tools_are_all_named(monkeypatch):
    monkeypatch.setattr(fetch.shutil, "which", lambda name: None)
    monkeypatch.setattr(fetch, "_eosvc", lambda: "")

    with pytest.raises(RuntimeError, match="conda, git, eosvc"):
        check_prerequisites()
