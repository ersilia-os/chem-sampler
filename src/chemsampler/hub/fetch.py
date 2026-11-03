"""Clone Hub models, download their checkpoints and build their conda environments.

Everything that starts a process (git, conda, eosvc, pip) goes through `_run`, the
one seam the tests replace.
"""

import datetime
import hashlib
import importlib.metadata
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

from ..store import state_path
from .envs import RDKIT_ENV, is_complete, read_marker, update_marker, write_marker

#: Models are cloned from here; the clone folder must be named exactly like the
#: model, because `eosvc` derives the S3 path from the folder name.
ORG_URL = "https://github.com/ersilia-os"

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
    #: Whether to build the `rdkit` environment too (`chemsampler fetch rdkit`),
    #: which the RDKit annotators run in.
    rdkit: bool = False
    #: The RDKit annotators the inputs list (`mw`, `qed`, `tanimoto`): they need
    #: that environment.
    rdkit_annotators: tuple[str, ...] = ()
    #: The version each model is fetched at: a release tag (e.g. "v1.0.0"), or
    #: `live` (its main branch's latest commit).
    versions: dict = field(default_factory=dict)


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
    """A clone's commit and whether it has local changes, after fetching its tags."""
    log.run(["git", "-C", model_dir, "fetch", "--quiet", "--tags"])
    changed = log.run(["git", "-C", model_dir, "status", "--porcelain"]).stdout
    commit = log.run(
        ["git", "-C", model_dir, "rev-parse", "--short", "HEAD"]
    ).stdout.strip()
    return RepoState(commit, bool(_local_changes(changed)))


#: The version that follows the repository's main branch, at its latest commit
#: (not necessarily a release).
LIVE = "live"

#: A release tag: v, then three numbers.
_RELEASE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")


def _release_key(tag: str) -> tuple[int, ...] | None:
    found = _RELEASE.match(tag)
    return None if found is None else tuple(int(n) for n in found.groups())


def _releases(lines) -> list[str]:
    """The release tags among `lines`, oldest first (compared as numbers)."""
    tags = {line.strip() for line in lines}
    return sorted((t for t in tags if _release_key(t)), key=_release_key)


def _commit_of(model_dir: Path, ref: str) -> str | None:
    found = _run(["git", "-C", model_dir, "rev-parse", "--verify", "--quiet", ref])
    return found.stdout.strip() if found.returncode == 0 else None


def _default_branch(model_dir: Path) -> str:
    """The remote's main branch, e.g. "main" (from origin/HEAD)."""
    found = _run(
        ["git", "-C", model_dir, "symbolic-ref", "--short", "refs/remotes/origin/HEAD"]
    )
    name = found.stdout.strip() if found.returncode == 0 else "origin/main"
    return name.split("/", 1)[-1]


def release_problem(model_dir: Path, tag: str) -> str | None:
    """
    Why a clone is not at version `tag`, or None if it is (or is no git clone,
    such as a model `ersilia fetch` left, which has nothing to compare). `live`
    asks only that the clone is on its main branch: whether that is the latest
    commit needs the network, which is `fetch`'s job.
    """
    if not (Path(model_dir) / ".git").is_dir():
        return None
    name = Path(model_dir).name
    if tag == LIVE:
        branch = _run(["git", "-C", model_dir, "symbolic-ref", "--short", "HEAD"])
        if branch.stdout.strip() != _default_branch(model_dir):
            return f"{name} is not at live (not on its main branch)"
        return None
    target = _commit_of(model_dir, f"refs/tags/{tag}^{{commit}}")
    if target is None:
        return f"{name}: its clone has no release {tag}"
    if _commit_of(model_dir, "HEAD") != target:
        return f"{name} is not at {tag}"
    return None


def _has_local_commits(model_dir: Path, log: _Log | None = None) -> bool:
    """Whether HEAD is on no remote branch and no tag: commits only this clone has.

    With no `log`, nothing is written down (what `chemsampler remove` wants)."""
    run = _run if log is None else log.run
    branches = run(["git", "-C", model_dir, "branch", "-r", "--contains", "HEAD"])
    tags = run(["git", "-C", model_dir, "tag", "--contains", "HEAD"])
    return not (branches.stdout.strip() or tags.stdout.strip())


def _to_live(
    model_id: str,
    model_dir: Path,
    state: RepoState,
    fresh: bool,
    log: _Log,
    result: "ModelResult",
    say: Callable[[str], None],
) -> bool:
    """
    Bring the clone to its main branch's latest commit, unless that would lose
    local work. Returns whether it is there afterwards.

    Raises
    ------
    BuildError
        If the checkout or the pull fails.
    """
    if fresh:  # a new clone is on the main branch, at its latest commit
        result.clone = f"cloned at {LIVE}"
        return True
    branch = _default_branch(model_dir)
    on_branch = release_problem(model_dir, LIVE) is None
    # commits only this clone has, on a branch or not, are never pulled over
    if state.dirty or _has_local_commits(model_dir, log):
        result.clone = f"not at {LIVE}"
        result.notes.append(f"not switched to {LIVE}: local changes or local commits")
        return False
    before = _commit_of(model_dir, "HEAD")
    if not on_branch:
        switched = log.run(["git", "-C", model_dir, "checkout", "--quiet", branch])
        if switched.returncode != 0:
            raise BuildError(f"git checkout {branch} failed: {_tail(switched)}")
    pulled = log.run(["git", "-C", model_dir, "pull", "--ff-only", "--quiet"])
    if pulled.returncode != 0:  # e.g. local commits on the branch itself
        result.clone = f"not at {LIVE}"
        result.notes.append(f"not updated to {LIVE}: {_tail(pulled)}")
        return False
    changed = _commit_of(model_dir, "HEAD") != before
    if changed:
        say(f"{model_id}: updated to {LIVE}")
    result.clone = f"updated to {LIVE}" if changed else f"at {LIVE}"
    return True


def _to_release(
    model_id: str,
    model_dir: Path,
    tag: str,
    state: RepoState,
    fresh: bool,
    log: _Log,
    result: "ModelResult",
    say: Callable[[str], None],
) -> bool:
    """
    Check the clone out at release `tag` unless that would lose local work.

    Returns whether it is at the tag afterwards.

    Raises
    ------
    BuildError
        If the repository has no such release, or the checkout fails.
    """
    target = _commit_of(model_dir, f"refs/tags/{tag}^{{commit}}")
    if target is None:
        tags = log.run(["git", "-C", model_dir, "tag", "-l"]).stdout.splitlines()
        known = ", ".join(_releases(tags)) or "none"
        raise BuildError(f"no release {tag} of {model_id}; its releases: {known}")
    if _commit_of(model_dir, "HEAD") == target:
        result.clone = f"cloned at {tag}" if fresh else f"at {tag}"
        return True
    if not fresh and (state.dirty or _has_local_commits(model_dir, log)):
        result.clone = f"not at {tag}"
        result.notes.append(f"not switched to {tag}: local changes or local commits")
        return False
    say(f"{model_id}: checking out {tag}")
    switched = log.run(["git", "-C", model_dir, "checkout", "--quiet", tag])
    if switched.returncode != 0:
        raise BuildError(f"git checkout {tag} failed: {_tail(switched)}")
    result.clone = f"cloned at {tag}" if fresh else f"switched to {tag}"
    return True


def _local_changes(porcelain: str) -> list[str]:
    """
    The lines of `git status --porcelain` that are the user's own work.

    Whatever eosvc leaves inside a clone is ours, not theirs, and must not make
    the clone look modified: that would stop fetch from ever switching its release.
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


# --- fetching -------------------------------------------------------------------------


def _load_state(models_dir: Path) -> dict:
    try:
        return json.loads(state_path(models_dir).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save_state(models_dir: Path, state: dict) -> None:
    state_path(models_dir).write_text(json.dumps(state, indent=2) + "\n")


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

    # 1. the clone, at the release the inputs name: cloned if absent, otherwise
    # checked out at that release (on request) unless it holds local work
    tag = plan.versions.get(model_id)
    if not tag:
        raise BuildError("no version (release tag) to fetch it at")
    fresh = not (model_dir / ".git").is_dir()
    if fresh:
        if model_dir.exists() and any(model_dir.iterdir()):
            raise BuildError(f"{model_dir} exists but is not a git clone")
        say(f"{model_id}: cloning")
        cloned = log.run(
            ["git", "clone", "--quiet", f"{ORG_URL}/{model_id}.git", model_dir]
        )
        if cloned.returncode != 0:
            raise BuildError(f"git clone failed: {_tail(cloned)}")
    state_now = repo_state(model_dir, log)
    if fresh or update:
        if tag == LIVE:
            _to_live(model_id, model_dir, state_now, fresh, log, result, say)
        else:
            _to_release(model_id, model_dir, tag, state_now, fresh, log, result, say)
        state_now = repo_state(model_dir, log)
    else:
        problem = release_problem(model_dir, tag)
        result.clone = f"at {tag}" if problem is None else f"not at {tag}"

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
    _ensure_env(result, env_dir, recipe, log, state_now.commit, rebuild, say)
    return result


def _ensure_env(
    result: ModelResult,
    env_dir: Path,
    recipe: Recipe,
    log: _Log,
    commit: str | None,
    rebuild: bool,
    say: Callable[[str], None],
) -> None:
    """Build an environment if absent, otherwise check it against its recipe."""
    name = result.model_id
    if env_dir.exists() and not is_complete(env_dir):
        say(f"{name}: removing an unfinished environment")
        shutil.rmtree(env_dir)
    if not env_dir.exists():
        say(f"{name}: building the environment (python {recipe.python})")
        build_env(env_dir, recipe, log, commit, name)
        result.env = "built"
        return
    check = check_env(env_dir, recipe, log)
    if check.matches:
        result.env = "ok"
        if commit and (read_marker(env_dir) or {}).get("commit") != commit:
            update_marker(env_dir, commit=commit)  # still right for this commit
    elif rebuild:
        say(f"{name}: the environment differs from its recipe, rebuilding")
        shutil.rmtree(env_dir)
        build_env(env_dir, recipe, log, commit, name)
        result.env = "rebuilt"
    else:
        result.env = "differs"
        result.differences = check.differences


def rdkit_recipe() -> Recipe:
    """
    The `rdkit` environment: this Python's version and the RDKit chemsampler runs
    with, so its scores equal those chemsampler computes itself.
    """
    python = f"{sys.version_info.major}.{sys.version_info.minor}"
    version = importlib.metadata.version("rdkit")
    text = f"python {python}\nrdkit=={version}\n"
    return Recipe(
        source="chemsampler",
        python=python,
        steps=(("pip", ["install", f"rdkit=={version}"]),),
        pip_pins={"rdkit": version},
        conda_pins={},
        unchecked=0,
        digest=hashlib.sha256(text.encode()).hexdigest(),
    )


def _fetch_rdkit(
    plan: FetchPlan, rebuild: bool, say: Callable[[str], None]
) -> ModelResult:
    """The `rdkit` environment, built or checked like a model's; it has no clone."""
    result = ModelResult(RDKIT_ENV, clone="not needed", checkpoints="none needed")
    log = _Log(plan.logs_dir / f"{RDKIT_ENV}.log")
    _ensure_env(
        result, plan.envs_dir / RDKIT_ENV, rdkit_recipe(), log, None, rebuild, say
    )
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

    Each model is cloned at the release tag `plan.versions` gives it. A failure
    on one model is recorded and does not stop the others. Existing clones are
    only inspected unless `update`, which checks them out at their release
    (never over local changes or local commits); existing environments are only
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
    if plan.rdkit:
        try:
            result = _fetch_rdkit(plan, rebuild, say)
        except (BuildError, OSError) as e:
            result = ModelResult(RDKIT_ENV, error=str(e))
            say(f"{RDKIT_ENV}: FAILED - {e}")
        results.append(result)
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
            "run `chemsampler fetch` to rebuild them."
        )
    unchecked = sum(r.unchecked for r in results if r.env in ("ok", "differs"))
    if unchecked:
        rows.append(
            f"{unchecked} recipe line(s) could not be compared (unpinned, GitHub or "
            "shell lines)."
        )
    return rows
