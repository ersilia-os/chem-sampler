"""Clone Hub models, download their checkpoints and build their conda environments.

Everything that starts a process (git, conda, eosvc, pip) goes through `_run`, the
one seam the tests replace.
"""

import datetime
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .envs import is_complete, write_marker

#: Models are cloned from here; the clone folder must be named exactly like the
#: model, because `eosvc` derives the S3 path from the folder name.
ORG_URL = "https://github.com/ersilia-os"

#: Remembers, next to the clones (not inside them, which would make them dirty),
#: for which commit the checkpoints were downloaded.
STATE_FILE = ".chemsampler-fetch.json"

#: What a model needs to run: its checkpoints. `fit` (the training data eosvc also
#: manages) is deliberately not downloaded.
CHECKPOINTS_PATH = "model/checkpoints"

#: eosvc marks a failure with this and exits non-zero even when a path it manages is
#: simply empty in S3 - which is the case for a model that keeps its data in git.
_EOSVC_FAILURE = "\u2717"
_EOSVC_EMPTY_PATH = "Nothing found at"

#: eosvc keeps its own state here, inside the clone. Model repositories do not
#: ignore it, so it would otherwise make every clone we fetched look modified.
_EOSVC_STATE_DIR = ".eosvc"

_LOW_DISK_GB = 10

#: `pip install` options that take a value, so the value is not read as a package.
_PIP_VALUE_OPTIONS = {
    "-i",
    "--index-url",
    "--extra-index-url",
    "-f",
    "--find-links",
    "-c",
    "--constraint",
    "-r",
    "--requirement",
    "-t",
    "--target",
    "--prefix",
    "--root",
    "--platform",
    "--python-version",
    "--implementation",
    "--abi",
    "--only-binary",
    "--no-binary",
}
_PY_TAG = re.compile(r"-py(\d)(\d+)")
_PIP_LINE = re.compile(r"^(?:python3?\s+-m\s+)?pip\s+(.*)$")


@dataclass(frozen=True)
class FetchPlan:
    """The models to fetch and where `chemsampler fetch` puts them."""

    model_ids: tuple[str, ...]
    models_dir: Path
    envs_dir: Path
    logs_dir: Path


@dataclass(frozen=True)
class Recipe:
    """How a model's environment is built, and what can be checked afterwards."""

    source: str  # "install.yml" or "Dockerfile"
    python: str  # "3.12"
    steps: tuple  # ("pip", args) | ("conda", args) | ("shell", command line)
    pip_pins: dict  # normalized package name -> version
    conda_pins: dict  # conda package name -> version
    unchecked: int  # lines whose outcome cannot be compared (unpinned, GitHub, shell)
    digest: str  # sha256 of the recipe file


@dataclass(frozen=True)
class EnvCheck:
    """Differences between an environment and its recipe."""

    differences: tuple[str, ...]

    @property
    def matches(self) -> bool:
        return not self.differences


@dataclass(frozen=True)
class RepoState:
    commit: str
    ahead: int
    behind: int
    dirty: bool


@dataclass
class ModelResult:
    """What happened to one model."""

    model_id: str
    clone: str = "-"
    checkpoints: str = "-"
    env: str = "-"
    differences: tuple[str, ...] = ()
    unchecked: int = 0
    error: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.error is None


class BuildError(RuntimeError):
    """A step of fetching a model failed."""


# --- processes -----------------------------------------------------------------------


def _run(argv, *, cwd=None, env=None) -> subprocess.CompletedProcess:
    """Run a command and capture its output. The only place that starts processes."""
    return subprocess.run(
        [str(a) for a in argv],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def _isolated_env() -> dict[str, str]:
    """Our environment with the user site off, so a build never leans on ~/.local."""
    return {**os.environ, "PYTHONNOUSERSITE": "1"}


class _Log:
    """Runs commands and appends each command and its output to one log file."""

    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)

    def run(self, argv, **kwargs) -> subprocess.CompletedProcess:
        result = _run(argv, **kwargs)
        with open(self.path, "a") as f:
            f.write(
                f"$ {shlex.join(str(a) for a in argv)}\n{result.stdout}{result.stderr}"
                f"[exit {result.returncode}]\n\n"
            )
        return result


def _tail(result: subprocess.CompletedProcess, size: int = 400) -> str:
    return (result.stderr or result.stdout or "").strip()[-size:]


def _eosvc_problem(output: str) -> str | None:
    """
    What went wrong in an eosvc run, or None if only empty paths were reported.

    eosvc exits non-zero when a path it manages holds nothing in S3, which is
    not a failure: the model simply keeps that data in git. Anything else it
    reports (no credentials, no network, a broken object) is a real failure.
    """
    failures = [
        line.strip()
        for line in output.splitlines()
        if line.strip().startswith(_EOSVC_FAILURE)
    ]
    if failures and all(_EOSVC_EMPTY_PATH in line for line in failures):
        return None
    return output.strip()[-400:] or "eosvc failed without saying why"


def _eosvc() -> str:
    sibling = Path(sys.executable).with_name("eosvc")
    return str(sibling) if sibling.exists() else (shutil.which("eosvc") or "")


def check_prerequisites() -> None:
    """Fail, naming everything missing, unless conda, git and eosvc can be run."""
    missing = [
        name
        for name, found in (
            ("conda", shutil.which("conda")),
            ("git", shutil.which("git")),
            ("eosvc", _eosvc()),
        )
        if not found
    ]
    if missing:
        raise RuntimeError(
            f"fetch needs {', '.join(missing)} on the PATH; install what is missing "
            "and run it again"
        )


def free_disk_gb(path: Path) -> float:
    """Free space on the disk that holds `path` (or its closest existing parent)."""
    probe = Path(path)
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe).free / 1e9


# --- recipes -------------------------------------------------------------------------


def _normalize(name: str) -> str:
    """The canonical form of a package name (PEP 503)."""
    return re.sub(r"[-_.]+", "-", name.split("[")[0]).lower()


def _version_matches(installed: str, pinned: str) -> bool:
    """
    Whether an installed version satisfies a recipe's pin, as pip judges it.

    A pin without a local label ignores the installed one (PEP 440), so
    `torch==2.9.1` is satisfied by the CPU wheel's `2.9.1+cpu` - which is what
    the models that install from the PyTorch CPU index actually get. A pin that
    does name a local label has to match it.
    """
    if "+" in pinned:
        return installed == pinned
    return installed.split("+", 1)[0] == pinned


def _is_url(spec: str) -> bool:
    return spec.startswith(("git+", "http://", "https://")) or "@" in spec


def _from_install_yml(path: Path, text: str) -> Recipe:
    spec = yaml.safe_load(text)
    steps, pip_pins, conda_pins, unchecked = [], {}, {}, 0
    for entry in spec["commands"]:
        tool, name, *rest = entry
        if tool == "pip":
            if _is_url(name) or not rest:
                steps.append(("pip", ["install", name, *rest]))
                unchecked += 1
            else:
                steps.append(("pip", ["install", f"{name}=={rest[0]}", *rest[1:]]))
                pip_pins[_normalize(name)] = str(rest[0])
        elif tool == "conda":
            rest = [x for x in rest if x not in ("-y", "--yes")]
            args = []
            if len(rest) >= 2:
                args += ["-c", str(rest[1])]
            if rest:
                args.append(f"{name}={rest[0]}")
                conda_pins[name] = str(rest[0])
            else:
                args.append(name)
                unchecked += 1
            steps.append(("conda", args))
        else:
            raise ValueError(f"install.yml: unknown command tool {tool!r}")
    return Recipe(
        source="install.yml",
        python=str(spec["python"]),
        steps=tuple(steps),
        pip_pins=pip_pins,
        conda_pins=conda_pins,
        unchecked=unchecked,
        digest=hashlib.sha256(text.encode()).hexdigest(),
    )


def _pip_packages(args: list[str]) -> list[str]:
    """The package specs in the arguments of a `pip install`, skipping options."""
    packages, skip = [], False
    for token in args[1:] if args and args[0] == "install" else args:
        if skip:
            skip = False
        elif token in _PIP_VALUE_OPTIONS:
            skip = True
        elif not token.startswith("-"):
            packages.append(token)
    return packages


def _from_dockerfile(path: Path, text: str) -> Recipe:
    from_line = next(
        (line for line in text.splitlines() if line.strip().startswith("FROM")), ""
    )
    match = _PY_TAG.search(from_line)
    if not match:
        raise ValueError(
            f"Dockerfile: could not find a python version in {from_line!r}"
        )
    steps, pip_pins, unchecked = [], {}, 0
    for line in text.replace("\\\n", " ").splitlines():
        line = line.strip()
        if not line.startswith("RUN "):
            continue
        command = line[len("RUN ") :].strip()
        pip = _PIP_LINE.match(command)
        try:
            args = shlex.split(pip.group(1)) if pip else None
        except ValueError:
            args = None
        if args and args[0] == "install":
            steps.append(("pip", args))
            for package in _pip_packages(args):
                if "==" in package and not _is_url(package):
                    name, version = package.split("==", 1)
                    pip_pins[_normalize(name)] = version.split(";")[0].strip()
                else:
                    unchecked += 1
        elif args:  # another pip command, e.g. upgrading pip itself
            steps.append(("pip", args))
        else:
            if re.match(r"^conda\s+install\b", command) and not re.search(
                r"(^|\s)(-y|--yes)(\s|$)", command
            ):
                command += " -y"
            steps.append(("shell", command))
            unchecked += 1
    return Recipe(
        source="Dockerfile",
        python=f"{match.group(1)}.{match.group(2)}",
        steps=tuple(steps),
        pip_pins=pip_pins,
        conda_pins={},
        unchecked=unchecked,
        digest=hashlib.sha256(text.encode()).hexdigest(),
    )


def read_recipe(model_dir) -> Recipe:
    """
    Read how a model's environment is built: `install.yml`, else the Dockerfile.

    Raises
    ------
    ValueError
        If the model has neither file, or the file cannot be understood.
    """
    model_dir = Path(model_dir)
    for name, parse in (
        ("install.yml", _from_install_yml),
        ("Dockerfile", _from_dockerfile),
    ):
        path = model_dir / name
        if path.is_file():
            try:
                return parse(path, path.read_text())
            except (KeyError, TypeError, ValueError, yaml.YAMLError) as e:
                raise ValueError(f"{path}: cannot read the recipe ({e})") from None
    raise ValueError(
        f"{model_dir}: no install.yml or Dockerfile to build an environment from"
    )


# --- environments --------------------------------------------------------------------


def build_env(
    env_dir: Path, recipe: Recipe, log: _Log, commit: str, model_id: str
) -> None:
    """
    Build a model's conda environment from its recipe, as declared.

    pip steps run through the environment's own python, never through whichever
    `pip` is first on PATH, which once installed into the wrong environment.
    The completion marker is written only after every step succeeded.

    Raises
    ------
    BuildError
        At the first step that fails; the unfinished environment is left behind
        without a marker, and the next fetch removes it.
    """
    env = _isolated_env()
    created = log.run(
        ["conda", "create", "-p", env_dir, f"python={recipe.python}", "-y"], env=env
    )
    if created.returncode != 0:
        raise BuildError(f"conda create failed: {_tail(created)}")
    python = env_dir / "bin" / "python"
    for number, (kind, payload) in enumerate(recipe.steps, start=1):
        if kind == "pip":
            argv = [python, "-m", "pip", *payload]
        elif kind == "conda":
            argv = ["conda", "install", "-p", env_dir, "-y", *payload]
        else:
            argv = ["conda", "run", "-p", env_dir, "bash", "-c", payload]
        result = log.run(argv, env=env)
        if result.returncode != 0:
            shown = payload if kind == "shell" else " ".join(payload)
            raise BuildError(
                f"step {number}/{len(recipe.steps)} ({kind}: {shown}) failed: {_tail(result)}"
            )
    write_marker(
        env_dir,
        model=model_id,
        commit=commit,
        built=datetime.datetime.now(tz=datetime.timezone.utc).isoformat(
            timespec="seconds"
        ),
        recipe=recipe.source,
        recipe_sha256=recipe.digest,
    )


def check_env(env_dir: Path, recipe: Recipe, log: _Log) -> EnvCheck:
    """
    Compare what is installed in a built environment with what its recipe pins.

    Checks the python version and every pinned pip and conda package, with one
    `pip list` (and one `conda list` only if the recipe has conda pins), comparing
    versions the way pip does (see `_version_matches`).
    """
    env = _isolated_env()
    python = env_dir / "bin" / "python"
    differences = []

    version = log.run(
        [python, "-c", "import sys; print('%d.%d' % sys.version_info[:2])"], env=env
    )
    if version.returncode != 0:
        return EnvCheck((f"the environment's python does not run: {_tail(version)}",))
    wanted = ".".join(recipe.python.split(".")[:2])
    if version.stdout.strip() != wanted:
        differences.append(
            f"python: installed {version.stdout.strip()}, recipe {wanted}"
        )

    listed = log.run(
        [python, "-m", "pip", "list", "--format=json", "--disable-pip-version-check"],
        env=env,
    )
    if listed.returncode != 0:
        return EnvCheck((*differences, f"pip list failed: {_tail(listed)}"))
    installed = {_normalize(p["name"]): p["version"] for p in json.loads(listed.stdout)}
    for name, pinned in recipe.pip_pins.items():
        if name not in installed:
            differences.append(f"{name}: not installed, recipe {pinned}")
        elif not _version_matches(installed[name], pinned):
            differences.append(f"{name}: installed {installed[name]}, recipe {pinned}")

    if recipe.conda_pins:
        listed = log.run(["conda", "list", "-p", env_dir, "--json"], env=env)
        if listed.returncode != 0:
            return EnvCheck((*differences, f"conda list failed: {_tail(listed)}"))
        conda_installed = {p["name"]: p["version"] for p in json.loads(listed.stdout)}
        for name, pinned in recipe.conda_pins.items():
            if name not in conda_installed:
                differences.append(f"{name} (conda): not installed, recipe {pinned}")
            elif not _version_matches(conda_installed[name], pinned):
                differences.append(
                    f"{name} (conda): installed {conda_installed[name]}, recipe {pinned}"
                )
    return EnvCheck(tuple(differences))


# --- clones ---------------------------------------------------------------------------


def repo_state(model_dir: Path, log: _Log) -> RepoState:
    """Where a clone stands against its remote, after a `git fetch` that changes no file."""
    log.run(["git", "-C", model_dir, "fetch", "--quiet"])
    counts = log.run(
        ["git", "-C", model_dir, "rev-list", "--left-right", "--count", "HEAD...@{u}"]
    )
    ahead = behind = 0
    if counts.returncode == 0 and len(counts.stdout.split()) == 2:
        ahead, behind = (int(n) for n in counts.stdout.split())
    changed = log.run(["git", "-C", model_dir, "status", "--porcelain"]).stdout
    dirty = bool(_local_changes(changed))
    commit = log.run(
        ["git", "-C", model_dir, "rev-parse", "--short", "HEAD"]
    ).stdout.strip()
    return RepoState(commit, ahead, behind, dirty)


def _local_changes(porcelain: str) -> list[str]:
    """
    The lines of `git status --porcelain` that are the user's own work.

    Whatever eosvc leaves inside a clone is ours, not theirs, and must not make
    the clone look modified: that would stop `--update` from ever pulling it.
    """
    lines = []
    for line in porcelain.splitlines():
        if not line.strip():
            continue
        path = line[3:].strip().strip('"')
        if path == _EOSVC_STATE_DIR or path.startswith(f"{_EOSVC_STATE_DIR}/"):
            continue
        lines.append(line)
    return lines


def _clone_note(state: RepoState) -> str:
    notes = []
    if state.behind and state.ahead:
        notes.append(f"diverged ({state.ahead} ahead, {state.behind} behind)")
    elif state.behind:
        notes.append(f"{state.behind} behind")
    elif state.ahead:
        notes.append(f"{state.ahead} ahead")
    if state.dirty:
        notes.append("local changes")
    return ", ".join(notes) or "up to date"


# --- fetching -------------------------------------------------------------------------


def _load_state(models_dir: Path) -> dict:
    try:
        return json.loads((models_dir / STATE_FILE).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save_state(models_dir: Path, state: dict) -> None:
    (models_dir / STATE_FILE).write_text(json.dumps(state, indent=2) + "\n")


def _fetch_one(
    model_id: str,
    plan: FetchPlan,
    state: dict,
    update: bool,
    rebuild: bool,
    say: Callable[[str], None],
) -> ModelResult:
    result = ModelResult(model_id)
    log = _Log(plan.logs_dir / f"{model_id}.log")
    model_dir = plan.models_dir / model_id
    env_dir = plan.envs_dir / model_id

    # 1. the clone: created if absent, otherwise only inspected (and updated on request)
    if not (model_dir / ".git").is_dir():
        if model_dir.exists() and any(model_dir.iterdir()):
            raise BuildError(f"{model_dir} exists but is not a git clone")
        say(f"{model_id}: cloning")
        cloned = log.run(
            ["git", "clone", "--quiet", f"{ORG_URL}/{model_id}.git", model_dir]
        )
        if cloned.returncode != 0:
            raise BuildError(f"git clone failed: {_tail(cloned)}")
        result.clone = "cloned"
        state_now = repo_state(model_dir, log)
    else:
        state_now = repo_state(model_dir, log)
        result.clone = _clone_note(state_now)
        if update:
            if state_now.dirty or state_now.ahead:
                result.notes.append("not updated: local changes or local commits")
            elif state_now.behind:
                say(f"{model_id}: updating {state_now.behind} commit(s)")
                pulled = log.run(
                    ["git", "-C", model_dir, "pull", "--ff-only", "--quiet"]
                )
                if pulled.returncode != 0:
                    raise BuildError(f"git pull failed: {_tail(pulled)}")
                state_now = repo_state(model_dir, log)
                result.clone = "updated"

    # 2. the checkpoints, which a clone alone does not carry
    if (model_dir / "access.json").is_file():
        recorded = state.get(model_id, {})
        if recorded.get("checkpoints_commit") == state_now.commit:
            result.checkpoints = "present"
        else:
            say(f"{model_id}: downloading checkpoints")
            downloaded = log.run(
                [_eosvc(), "download", "--path", CHECKPOINTS_PATH], cwd=model_dir
            )
            if downloaded.returncode == 0:
                result.checkpoints = "downloaded"
            else:
                problem = _eosvc_problem(
                    (downloaded.stdout or "") + (downloaded.stderr or "")
                )
                if problem is not None:
                    raise BuildError(f"eosvc download failed: {problem}")
                result.checkpoints = "none in S3"
                result.notes.append(
                    f"nothing at {CHECKPOINTS_PATH} in S3; its data must be in git"
                )
            state[model_id] = {**recorded, "checkpoints_commit": state_now.commit}
    else:
        result.checkpoints = "none needed"

    # 3. the environment: built if absent, otherwise checked against the recipe
    recipe = read_recipe(model_dir)
    result.unchecked = recipe.unchecked
    if env_dir.exists() and not is_complete(env_dir):
        say(f"{model_id}: removing an unfinished environment")
        shutil.rmtree(env_dir)
    if not env_dir.exists():
        say(f"{model_id}: building the environment (python {recipe.python})")
        build_env(env_dir, recipe, log, state_now.commit, model_id)
        result.env = "built"
        return result
    check = check_env(env_dir, recipe, log)
    if check.matches:
        result.env = "ok"
    elif rebuild:
        say(f"{model_id}: the environment differs from its recipe, rebuilding")
        shutil.rmtree(env_dir)
        build_env(env_dir, recipe, log, state_now.commit, model_id)
        result.env = "rebuilt"
    else:
        result.env = "differs"
        result.differences = check.differences
    return result


def fetch_models(
    plan: FetchPlan,
    *,
    update: bool = False,
    rebuild: bool = False,
    say: Callable[[str], None] = print,
) -> list[ModelResult]:
    """
    Clone, download and build every model of the plan, one at a time.

    A failure on one model is recorded and does not stop the others. Existing
    clones are only inspected unless `update`; existing environments are only
    checked against their recipe unless `rebuild`, which rebuilds the ones that
    differ.
    """
    for folder in (plan.models_dir, plan.envs_dir, plan.logs_dir):
        folder.mkdir(parents=True, exist_ok=True)
    free = free_disk_gb(plan.envs_dir)
    if free < _LOW_DISK_GB:
        say(f"warning: only {free:.1f} GB free where the environments go")
    state = _load_state(plan.models_dir)
    results = []
    for model_id in plan.model_ids:
        try:
            result = _fetch_one(model_id, plan, state, update, rebuild, say)
        except (BuildError, ValueError, OSError) as e:
            result = ModelResult(model_id, error=str(e))
            say(f"{model_id}: FAILED - {e}")
        results.append(result)
        _save_state(plan.models_dir, state)
    return results


def summarize(results: list[ModelResult]) -> list[str]:
    """The report `fetch` prints at the end: a table, then what needs attention."""
    width = max([len("model"), *(len(r.model_id) for r in results)])
    rows = [f"{'model':<{width}}  {'clone':<22} {'checkpoints':<13} environment"]
    for r in results:
        rows.append(
            f"{r.model_id:<{width}}  {r.clone:<22} {r.checkpoints:<13} "
            f"{'FAILED' if r.error else r.env}"
        )
    for r in results:
        for note in r.notes:
            rows.append(f"{r.model_id}: {note}")
        if r.error:
            rows.append(f"{r.model_id}: {r.error}")
        for difference in r.differences:
            rows.append(f"{r.model_id}: environment differs - {difference}")
    differing = [r.model_id for r in results if r.env == "differs"]
    if differing:
        rows.append(
            f"{len(differing)} environment(s) differ from their recipe ({', '.join(differing)}); "
            "run `chemsampler fetch --rebuild` to rebuild them."
        )
    unchecked = sum(r.unchecked for r in results if r.env in ("ok", "differs"))
    if unchecked:
        rows.append(
            f"{unchecked} recipe line(s) could not be compared (unpinned, GitHub or "
            "shell lines)."
        )
    return rows
