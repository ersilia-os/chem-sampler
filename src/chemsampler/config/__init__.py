"""Load and validate chemsampler's inputs: similarity, generators and annotators."""

import csv
import datetime
from dataclasses import dataclass
from pathlib import Path

from .. import store
from ..hub.catalog import validate_models
from ..hub.client import Backend
from ..hub.envs import env_path, is_complete
from ..hub.fetch import FetchPlan
from ..models.annotator import RDKIT_ANNOTATORS
from ..models.generator import GeneratorPool, HubGenerator
from ..models.hub_annotator import HubAnnotator
from ..models.similarity import HubSimilarity
from ..models.spec import AnnotatorSpec
from ..project import PROJECT_FILE, Project, find_project
from ..utils.jsonfile import is_list_of, is_text, read_object
from ..utils.logging import logger
from .schema import (
    OPTIONAL_FIELDS,
    REQUIRED_FIELDS,
    SOURCES,
    AnnotatorRow,
    check_header,
    parse_annotator,
)

#: The three model lists. Each comes from its flag (a CSV file), a key of the config
#: file, or a CSV of the same name in the `path_to_models` folder (or, failing
#: all of those, in the working directory). The CSVs' single column has the same
#: name; annotators.csv has the columns of `schema.REQUIRED_FIELDS`.
INPUT_NAMES = ("similarity", "generators", "annotators")
PATH_KEYS = (
    "path_to_models",
    "path_to_model_repositories",
    "path_to_envs",
    "path_to_rounds",
)
CONFIG_KEYS = (*PATH_KEYS, *INPUT_NAMES)


@dataclass(frozen=True)
class Inputs:
    """Everything `hill_climb` needs from the user, validated and built."""

    generators: GeneratorPool | None
    similarity: HubSimilarity | None
    annotators: list[AnnotatorSpec]
    models_dir: str | None
    envs_dir: str | None
    rounds_dir: str | None

    def describe(self) -> dict:
        """What was resolved, to save next to a run's results."""
        return {
            "similarity": None if self.similarity is None else self.similarity.model_id,
            "generators": (
                None
                if self.generators is None
                else [g.model_id for g in self.generators.generators]
            ),
            "annotators": [
                {
                    "model": spec.annotator_id,
                    "cutoff": list(spec.cutoff)
                    if isinstance(spec.cutoff, tuple)
                    else spec.cutoff,
                    "direction": spec.direction,
                    "priority": spec.priority,
                    "weight": spec.weight,
                }
                for spec in self.annotators
            ],
            "path_to_model_repositories": self.models_dir,
            "path_to_envs": self.envs_dir,
        }

    def new_round_dir(self) -> Path | None:
        """A fresh, dated folder under `rounds/` for one run's results."""
        if self.rounds_dir is None:
            return None
        # Local time: the folder name is read by whoever started the run.
        stamp = datetime.datetime.now().astimezone().strftime("%Y-%m-%d_%H-%M-%S")
        folder = Path(self.rounds_dir) / stamp
        suffix = 2
        while folder.exists():
            folder = Path(self.rounds_dir) / f"{stamp}_{suffix}"
            suffix += 1
        return folder


def load_inputs(
    config: str | None = None,
    path_to_models: str | None = None,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
    *,
    similarity: str | None = None,
    generators: str | None = None,
    annotators: str | None = None,
    path_to_envs: str | None = None,
) -> Inputs:
    """
    Load the inputs and check them against the catalog of Ersilia Hub models.

    Annotators are required, and so is at least one of similarity and
    generators. Each input comes from the first of these that has it: its own
    CSV file (`similarity`, `generators`, `annotators`); the config file's key of
    that name, or `<path_to_models>/<name>.csv` (both is an error); `<name>.csv`
    in the working directory, which is looked at only when no folder is named.

    Parameters
    ----------
    config : str, optional
        JSON file with the keys `path_to_models`, `path_to_model_repositories`,
        `similarity`, `generators` and `annotators`. Relative paths in it are
        resolved against the file's folder.
    path_to_models : str, optional
        Folder with `similarity.csv`, `generators.csv` and `annotators.csv`.
        Overrides the config file's key.
    path_to_model_repositories : str, optional
        Folder with one sub-folder per model. Required by the "run_sh" backend
        and overrides the config file's key.
    backend : {"ersilia", "run_sh"}, optional
        Passed to every Hub model, by default "run_sh".
    similarity, generators, annotators : str, optional
        CSV files with that input, each with its own header (see the README).
    path_to_envs : str, optional
        Folder with one conda environment per model, for models cloned and built
        by `chemsampler fetch`. Overrides the config file's key. When it is set,
        every model needs a complete environment there.

    Notes
    -----
    Inside a project (a folder with `chemsampler.json`, here or in a parent) and
    with no `config`, the project's marker file is the config. Models are looked
    for in the machine's store unless the flags, the config or the marker name a
    folder.

    Returns
    -------
    Inputs
        The generator pool and the similarity model (each `None` if not given),
        the annotators and the folder of model repositories.

    Raises
    ------
    ValueError
        If annotators, or both similarity and generators, are missing; an input
        is given twice or malformed; a model is not in the catalog with the label
        its role needs; or a model folder is missing.
    """
    files = {
        "similarity": similarity,
        "generators": generators,
        "annotators": annotators,
    }
    loaded = _load(
        INPUT_NAMES,
        config,
        path_to_models,
        path_to_model_repositories,
        backend,
        files,
        path_to_envs,
    )
    return Inputs(
        generators=loaded["generators"],
        similarity=loaded["similarity"],
        annotators=loaded["annotators"],
        models_dir=loaded["models_dir"],
        envs_dir=loaded["envs_dir"],
        rounds_dir=loaded["rounds_dir"],
    )


def load_annotators(
    config: str | None = None,
    path_to_models: str | None = None,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
    *,
    annotators: str | None = None,
    path_to_envs: str | None = None,
) -> list[AnnotatorSpec]:
    """
    Load only the annotators, for commands that score without generating.

    Takes the same arguments as `load_inputs` and applies the same rules, but
    reads and checks nothing else. The repositories folder is only needed when an
    annotator is an Ersilia model.

    Returns
    -------
    list[AnnotatorSpec]
        One spec per annotator, in file order.
    """
    loaded = _load(
        ("annotators",),
        config,
        path_to_models,
        path_to_model_repositories,
        backend,
        {"annotators": annotators},
        path_to_envs,
    )
    return loaded["annotators"]


def load_annotator(
    source: str,
    model: str,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
    *,
    path_to_envs: str | None = None,
):
    """
    Build and check one annotator, to score a molecule against a single model.

    An Ersilia model gets the same checks as in `load_annotators` (repositories
    folder for "run_sh", catalog label, model folder and environment); an rdkit
    one needs none.

    Parameters
    ----------
    source : {"ersilia", "rdkit"}
        Where the annotator comes from.
    model : str
        An Ersilia Hub model id, or an in-house name ("mw", "qed").
    path_to_model_repositories, path_to_envs : str, optional
        Folders with the model repositories and their environments. They
        default to the machine's store.
    backend : {"ersilia", "run_sh"}, optional
        By default "run_sh".

    Returns
    -------
    MolecularWeightAnnotator, QEDAnnotator or HubAnnotator
        An object with a `score(smiles_list)` method.
    """
    models_dir = envs_dir = None
    if source == "ersilia":
        loc = _locations(None, None, path_to_model_repositories, path_to_envs)
        envs_dir = _as_text(loc.envs)
        if backend == "run_sh":
            models_dir = _check_repositories(loc.repos, loc.from_store)
        else:
            models_dir = _as_text(loc.repos)
        validate_models(annotators=[model])
        if backend == "run_sh":
            _check_model_folders(models_dir, envs_dir, [model])
    return build_annotator(source, model, None, backend, models_dir, envs_dir)


def build_annotator(
    source: str,
    model: str,
    column: str | None = None,
    backend: Backend = "run_sh",
    models_dir: str | None = None,
    envs_dir: str | None = None,
):
    """
    Build an annotator from its source and model name.

    Parameters
    ----------
    source : {"ersilia", "rdkit"}
        Where the annotator comes from.
    model : str
        An Ersilia Hub model id for "ersilia", or an in-house name for "rdkit"
        ("mw", "qed").
    column : str, optional
        Output column to score on, for Ersilia models with several numeric
        outputs. Not used by rdkit annotators.
    backend : {"ersilia", "run_sh"}, optional
        Passed to `HubAnnotator`, by default "run_sh". Not used by rdkit.
    models_dir : str, optional
        Folder with the model repositories, passed to `HubAnnotator`.
    envs_dir : str, optional
        Folder with the models' conda environments, passed to `HubAnnotator`.

    Returns
    -------
    MolecularWeightAnnotator, QEDAnnotator or HubAnnotator
        An object with a `score(smiles_list)` method.

    Raises
    ------
    ValueError
        If `source` is unknown, or `model` is not an in-house annotator when
        `source` is "rdkit".
    """
    if source == "rdkit":
        if model not in RDKIT_ANNOTATORS:
            raise ValueError(
                f"unknown rdkit model {model!r}; choose from {sorted(RDKIT_ANNOTATORS)}"
            )
        return RDKIT_ANNOTATORS[model]()
    if source == "ersilia":
        return HubAnnotator(
            model,
            column=column,
            backend=backend,
            models_dir=models_dir,
            envs_dir=envs_dir,
        )
    raise ValueError(f"source must be one of {list(SOURCES)}, got {source!r}")


# --- reading ----------------------------------------------------------------------


def _read_config(path: Path) -> dict:
    """Read and type-check a config file."""
    data = read_object(path, expected=f" with keys {list(CONFIG_KEYS)}")
    unknown = sorted(set(data) - set(CONFIG_KEYS))
    if unknown:
        raise ValueError(
            f"{path}: unknown key(s) {unknown}; allowed: {list(CONFIG_KEYS)}"
        )
    for key in PATH_KEYS:
        if key in data and not is_text(data[key]):
            raise ValueError(f"{path}: '{key}' must be a path (text)")
    for key in ("similarity", "generators"):
        if key in data and not is_list_of(data[key], str):
            raise ValueError(f"{path}: '{key}' must be a list of model ids")
    if "annotators" in data and not is_list_of(data["annotators"], dict):
        raise ValueError(f"{path}: 'annotators' must be a list of objects")
    return data


def _unique_ids(ids: list[str], where: str) -> list[str]:
    """Strip the ids and require at least one, with no blanks or repeats."""
    ids = [i.strip() for i in ids]
    if not ids:
        raise ValueError(f"{where}: at least one model is required")
    if "" in ids:
        raise ValueError(f"{where}: empty model id")
    repeated = sorted({i for i in ids if ids.count(i) > 1})
    if repeated:
        raise ValueError(f"{where}: listed more than once: {repeated}")
    return ids


def _read_ids_csv(path: Path, column: str) -> list[str]:
    """A one-column CSV whose header is the input's name; blank rows are skipped."""
    with open(path, newline="") as f:
        reader = csv.reader(f)
        check_header(path, next(reader, None), [column])
        ids = []
        for number, row in enumerate(reader, start=1):
            if not any(cell.strip() for cell in row):
                continue
            if len(row) != 1:
                raise ValueError(
                    f"{path.name} row {number}: expected one value, got {len(row)}"
                )
            ids.append(row[0])
    return _unique_ids(ids, path.name)


def _read_annotators_csv(path: Path) -> list[AnnotatorRow]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        check_header(path, reader.fieldnames, REQUIRED_FIELDS, OPTIONAL_FIELDS)
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        rows = []
        for number, raw in enumerate(reader, start=1):
            if raw.get(None):
                raise ValueError(
                    f"{path.name} row {number}: more values than header columns"
                )
            row = {k: v for k, v in raw.items() if k is not None}
            if not any((v or "").strip() for v in row.values()):
                continue
            rows.append(parse_annotator(row, f"{path.name} row {number}"))
    return rows


def _find_source(name, flag, data, models_path, read_csv, read_inline, required):
    """
    The contents of one input, from the first place that has it.

    In order: its flag; the config file's inline key, or its CSV in the
    `path_to_models` folder (both is an error); its CSV in the working directory,
    only when no folder was named. Returns None for an optional input nobody has.
    """
    if flag is not None:
        return read_csv(Path(flag).expanduser())
    csv_path = models_path / f"{name}.csv" if models_path is not None else None
    has_csv = csv_path is not None and csv_path.is_file()
    if name in data and has_csv:
        raise ValueError(
            f"{name} is given twice: in the config file and in {csv_path}; use one"
        )
    if name in data:
        return read_inline(data[name])
    if has_csv:
        return read_csv(csv_path)
    if models_path is None:
        in_cwd = Path.cwd() / f"{name}.csv"
        if in_cwd.is_file():
            logger.info(f"Using {name}.csv from the working directory")
            return read_csv(in_cwd)
    if not required:
        return None
    where = f" (looked in {models_path})" if models_path is not None else ""
    raise ValueError(
        f"{name} is required: pass --{name} FILE, add the key '{name}' to the "
        f"config file, put {name}.csv in --path-to-models{where}, or have "
        f"{name}.csv in the working directory"
    )


def _annotator_rows(items: list[dict]) -> list[AnnotatorRow]:
    if not items:
        raise ValueError("annotators: at least one annotator is required")
    return [
        parse_annotator(item, f"config file annotators[{i}]")
        for i, item in enumerate(items)
    ]


#: The `Project` attribute that holds each folder, for configs that are project markers.
_PROJECT_FOLDER = {
    "path_to_models": "input",
    "path_to_model_repositories": "models",
    "path_to_envs": "envs",
    "path_to_rounds": "rounds",
}


@dataclass(frozen=True)
class _Locations:
    """The config contents and the folders, flags taking precedence over the config."""

    data: dict
    models: Path | None
    repos: Path | None
    envs: Path | None
    rounds: Path | None
    root: Path | None  # the project root, when the config is a project marker
    from_store: (
        bool  # the model repositories are the store's, not a folder the user named
    )


def _locations(
    config, path_to_models, path_to_model_repositories, path_to_envs=None
) -> _Locations:
    """
    Read the config and settle the folders.

    With no config given, the marker file of the project the working directory is
    in (if any) is the config, and the folders it does not name get their
    defaults. Relative paths in a config resolve against the config file.

    The model repositories and the environments come from the first of: a flag,
    the config, or the machine's store. The store is a pair: its environments are
    used only together with its model repositories, so a repositories folder the
    user names alone is taken to hold ready-to-run models (as `ersilia fetch`
    leaves them), without environments.
    """
    root = project = None
    if config is None:
        root = find_project()
        if root is not None:
            config = root / PROJECT_FILE
            project = Project.load(root)
    data, base = {}, Path.cwd()
    if config is not None:
        config_path = Path(config)
        data = _read_config(config_path)
        base = config_path.resolve().parent

    def folder(flag, key):
        if flag is not None:
            return Path(flag).expanduser()
        if key in data:
            path = Path(data[key]).expanduser()
            return path if path.is_absolute() else base / path
        if project is not None:
            return getattr(project, _PROJECT_FOLDER[key])
        return None

    repos = folder(path_to_model_repositories, "path_to_model_repositories")
    envs = folder(path_to_envs, "path_to_envs")
    from_store = repos is None
    if from_store:
        shared = store.resolve()
        repos = shared.models
        envs = envs if envs is not None else shared.envs
    return _Locations(
        data,
        folder(path_to_models, "path_to_models"),
        repos,
        envs,
        folder(None, "path_to_rounds"),
        root,
        from_store,
    )


# --- loading ----------------------------------------------------------------------


def _read_inputs(names, loc: _Locations, files):
    """Read the named inputs and check them without looking at the Hub or the disk."""
    if loc.models is not None and not loc.models.is_dir():
        raise ValueError(f"--path-to-models {loc.models} is not a folder")

    def ids(name):
        return _find_source(
            name,
            files.get(name),
            loc.data,
            loc.models,
            lambda p: _read_ids_csv(p, name),
            lambda values: _unique_ids(values, name),
            required=False,
        )

    generators = similarity = None
    if "generators" in names:
        generators = ids("generators")
    if "similarity" in names:
        similarity = ids("similarity")
        if similarity is not None and len(similarity) != 1:
            raise ValueError(
                "similarity: only one similarity model is supported for now, "
                f"got {len(similarity)}"
            )
    if "generators" in names and generators is None and similarity is None:
        raise ValueError(
            "give at least one of similarity and generators: pass --similarity or "
            "--generators FILE, add the key to the config file, put similarity.csv "
            "or generators.csv in --path-to-models, or have one in the working "
            "directory"
        )
    rows = _find_source(
        "annotators",
        files.get("annotators"),
        loc.data,
        loc.models,
        _read_annotators_csv,
        _annotator_rows,
        required=True,
    )
    _check_annotator_rows(rows)
    return generators, similarity, rows


def _load(
    names,
    config,
    path_to_models,
    path_to_model_repositories,
    backend,
    files,
    path_to_envs=None,
) -> dict:
    """
    Read, validate and build the named inputs.

    The order is cheap checks first: file contents, then the repositories
    folder, then the catalog, then each model's own folder and environment.
    """
    loc = _locations(config, path_to_models, path_to_model_repositories, path_to_envs)
    generators, similarity, rows = _read_inputs(names, loc, files)

    ersilia_annotators = [row.model for row in rows if row.source == "ersilia"]
    hub_ids = [*(generators or []), *(similarity or []), *ersilia_annotators]

    run_sh = backend == "run_sh" and bool(hub_ids)
    models_dir = (
        _check_repositories(loc.repos, loc.from_store)
        if run_sh
        else _as_text(loc.repos)
    )
    envs_dir = _as_text(loc.envs)

    validate_models(
        generators=generators or [],
        similarity=similarity or [],
        annotators=ersilia_annotators,
    )
    if run_sh:
        _check_model_folders(models_dir, envs_dir, hub_ids)

    loaded = {
        "models_dir": models_dir,
        "envs_dir": envs_dir,
        "rounds_dir": _as_text(loc.rounds),
        "annotators": [_to_spec(row, backend, models_dir, envs_dir) for row in rows],
        "generators": None,
        "similarity": None,
    }
    if generators is not None:
        loaded["generators"] = GeneratorPool(
            [
                HubGenerator(
                    i, backend=backend, models_dir=models_dir, envs_dir=envs_dir
                )
                for i in generators
            ]
        )
    if similarity is not None:
        loaded["similarity"] = HubSimilarity(
            similarity[0], backend=backend, models_dir=models_dir, envs_dir=envs_dir
        )
    return loaded


def locate_models(
    config: str | None = None,
    path_to_model_repositories: str | None = None,
    path_to_envs: str | None = None,
) -> tuple[Path, Path | None]:
    """
    Where the models and their environments are, as `run` and `fetch` see them.

    The flags, then the config (or the project's marker file), then the machine's
    store. The environments folder is None for a repositories folder the user
    named alone, which holds models that bring their own interpreter.

    Returns
    -------
    tuple[Path, Path or None]
        The folder with the clones and the one with the environments.
    """
    loc = _locations(config, None, path_to_model_repositories, path_to_envs)
    return loc.repos, loc.envs


def plan_fetch(
    config: str | None = None,
    path_to_models: str | None = None,
    path_to_model_repositories: str | None = None,
    path_to_envs: str | None = None,
    *,
    similarity: str | None = None,
    generators: str | None = None,
    annotators: str | None = None,
) -> FetchPlan:
    """
    List the Hub models the inputs need and the folders to put them in.

    Reads and checks the inputs like `load_inputs` does, but looks at no model
    folder, since creating them is what `chemsampler fetch` is for. Models that
    appear in several inputs are listed once.

    The clones go to the model repositories folder and the environments to the
    envs folder: the flags, the config, or else the machine's store. Logs go in
    the store's `logs` folder, or in a `logs` folder beside the clones when the
    user named the folder for them.

    Raises
    ------
    ValueError
        If the inputs are invalid.
    """
    loc = _locations(config, path_to_models, path_to_model_repositories, path_to_envs)
    files = {
        "similarity": similarity,
        "generators": generators,
        "annotators": annotators,
    }
    generators_ids, similarity_ids, rows = _read_inputs(INPUT_NAMES, loc, files)
    ersilia_annotators = [row.model for row in rows if row.source == "ersilia"]
    validate_models(
        generators=generators_ids or [],
        similarity=similarity_ids or [],
        annotators=ersilia_annotators,
    )
    shared = store.resolve()
    envs_dir = loc.envs if loc.envs is not None else shared.envs
    wanted = [*(generators_ids or []), *(similarity_ids or []), *ersilia_annotators]
    return FetchPlan(
        model_ids=tuple(dict.fromkeys(wanted)),
        models_dir=loc.repos,
        envs_dir=envs_dir,
        logs_dir=shared.logs if loc.from_store else loc.repos.parent / "logs",
    )


def _check_annotator_rows(rows: list[AnnotatorRow]) -> None:
    models = [row.model for row in rows]
    repeated = sorted({m for m in models if models.count(m) > 1})
    if repeated:
        raise ValueError(f"annotators: listed more than once: {repeated}")
    if not any(row.priority == "primary" for row in rows):
        raise ValueError(
            "annotators: at least one primary annotator is required; secondary "
            "annotators are guard-rails and are not optimized"
        )


def _as_text(path: Path | None) -> str | None:
    return None if path is None else str(path)


def _check_repositories(repos: Path, from_store: bool) -> str:
    """The folder with the model repositories, which must exist."""
    if not repos.is_dir():
        if from_store:
            raise ValueError(
                f"no models fetched yet (looked in {repos}); run `chemsampler fetch`"
            )
        raise ValueError(f"the model repositories folder {repos} does not exist")
    return str(repos)


def _check_model_folders(
    models_dir: str, envs_dir: str | None, model_ids: list[str]
) -> None:
    """Every model needs its folder and, when environments are used, a built one."""
    missing = [i for i in model_ids if not (Path(models_dir) / i).is_dir()]
    if missing:
        hint = (
            "run `chemsampler fetch` to clone them"
            if envs_dir
            else "the run_sh backend needs one folder per model there"
        )
        raise ValueError(f"no folder for {missing} in {models_dir}; {hint}")
    if envs_dir:
        unbuilt = [i for i in model_ids if not is_complete(env_path(envs_dir, i))]
        if unbuilt:
            raise ValueError(
                f"no complete environment for {unbuilt} in {envs_dir}; run "
                "`chemsampler fetch` to build them"
            )


def _to_spec(
    row: AnnotatorRow, backend: Backend, models_dir: str | None, envs_dir: str | None
):
    return AnnotatorSpec(
        annotator_id=row.model,
        annotator=build_annotator(
            row.source, row.model, row.column, backend, models_dir, envs_dir
        ),
        cutoff=row.cutoff,
        direction=row.direction,
        weight=1.0 if row.weight is None else row.weight,
        priority=row.priority,
    )
