"""Load and validate chemsampler's inputs: similarity, generators and annotators."""

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from ..hub.catalog import validate_models
from ..hub.client import Backend
from ..models.annotator import RDKIT_ANNOTATORS
from ..models.generator import GeneratorPool, HubGenerator
from ..models.hub_annotator import HubAnnotator
from ..models.similarity import HubSimilarity
from ..models.spec import AnnotatorSpec
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
CONFIG_KEYS = ("path_to_models", "path_to_model_repositories", *INPUT_NAMES)


@dataclass(frozen=True)
class Inputs:
    """Everything `hill_climb` needs from the user, validated and built."""

    generators: GeneratorPool | None
    similarity: HubSimilarity | None
    annotators: list[AnnotatorSpec]
    models_dir: str | None


def load_inputs(
    config: str | None = None,
    path_to_models: str | None = None,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
    *,
    similarity: str | None = None,
    generators: str | None = None,
    annotators: str | None = None,
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
    )
    return Inputs(
        generators=loaded["generators"],
        similarity=loaded["similarity"],
        annotators=loaded["annotators"],
        models_dir=loaded["models_dir"],
    )


def load_annotators(
    config: str | None = None,
    path_to_models: str | None = None,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
    *,
    annotators: str | None = None,
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
    )
    return loaded["annotators"]


def load_annotator(
    source: str,
    model: str,
    path_to_model_repositories: str | None = None,
    backend: Backend = "run_sh",
):
    """
    Build and check one annotator, to score a molecule against a single model.

    An Ersilia model gets the same checks as in `load_annotators` (repositories
    folder for "run_sh", Hub label, model folder); an rdkit one needs none.

    Parameters
    ----------
    source : {"ersilia", "rdkit"}
        Where the annotator comes from.
    model : str
        An Ersilia Hub model id, or an in-house name ("mw", "qed").
    path_to_model_repositories : str, optional
        Folder with the model repositories; required for an Ersilia model with
        the "run_sh" backend.
    backend : {"ersilia", "run_sh"}, optional
        By default "run_sh".

    Returns
    -------
    MolecularWeightAnnotator, QEDAnnotator or HubAnnotator
        An object with a `score(smiles_list)` method.
    """
    models_dir = None
    if source == "ersilia":
        repos = path_to_model_repositories
        if backend == "run_sh":
            models_dir = _check_repositories(None if repos is None else Path(repos))
        else:
            models_dir = repos
        validate_models(annotators=[model])
        if backend == "run_sh":
            _check_model_folders(models_dir, [model])
    return build_annotator(source, model, None, backend, models_dir)


def build_annotator(
    source: str,
    model: str,
    column: str | None = None,
    backend: Backend = "run_sh",
    models_dir: str | None = None,
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
            model, column=column, backend=backend, models_dir=models_dir
        )
    raise ValueError(f"source must be one of {list(SOURCES)}, got {source!r}")


# --- reading ----------------------------------------------------------------------


def _read_config(path: Path) -> dict:
    """Read and type-check a config file."""
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: not valid JSON ({e})") from None
    if not _is_object(data):
        raise ValueError(
            f"{path}: expected a JSON object with keys {list(CONFIG_KEYS)}"
        )
    unknown = sorted(set(data) - set(CONFIG_KEYS))
    if unknown:
        raise ValueError(
            f"{path}: unknown key(s) {unknown}; allowed: {list(CONFIG_KEYS)}"
        )
    for key in ("path_to_models", "path_to_model_repositories"):
        if key in data and not isinstance(data[key], str):
            raise ValueError(f"{path}: '{key}' must be a path (text)")
    for key in ("similarity", "generators"):
        if key in data and not _is_list_of(data[key], str):
            raise ValueError(f"{path}: '{key}' must be a list of model ids")
    if "annotators" in data and not _is_list_of(data["annotators"], dict):
        raise ValueError(f"{path}: 'annotators' must be a list of objects")
    return data


def _is_object(value) -> bool:
    return isinstance(value, dict)


def _is_list_of(value, kind) -> bool:
    return isinstance(value, list) and all(isinstance(item, kind) for item in value)


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


def _locations(config, path_to_models, path_to_model_repositories):
    """The config contents and the two folders, flags taking precedence."""
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
        return None

    return (
        data,
        folder(path_to_models, "path_to_models"),
        folder(path_to_model_repositories, "path_to_model_repositories"),
    )


# --- loading ----------------------------------------------------------------------


def _load(
    names, config, path_to_models, path_to_model_repositories, backend, files
) -> dict:
    """
    Resolve, validate and build the named inputs.

    The order is cheap checks first: file contents, then the repositories
    folder, then the catalog, then each model's own folder.
    """
    data, models_path, repos = _locations(
        config, path_to_models, path_to_model_repositories
    )
    if models_path is not None and not models_path.is_dir():
        raise ValueError(f"--path-to-models {models_path} is not a folder")

    def ids(name):
        return _find_source(
            name,
            files.get(name),
            data,
            models_path,
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
        data,
        models_path,
        _read_annotators_csv,
        _annotator_rows,
        required=True,
    )
    _check_annotator_rows(rows)

    ersilia_annotators = [row.model for row in rows if row.source == "ersilia"]
    hub_ids = [*(generators or []), *(similarity or []), *ersilia_annotators]

    run_sh = backend == "run_sh" and bool(hub_ids)
    models_dir = _check_repositories(repos) if run_sh else _as_text(repos)

    validate_models(
        generators=generators or [],
        similarity=similarity or [],
        annotators=ersilia_annotators,
    )
    if run_sh:
        _check_model_folders(models_dir, hub_ids)

    loaded = {
        "models_dir": models_dir,
        "annotators": [_to_spec(row, backend, models_dir) for row in rows],
        "generators": None,
        "similarity": None,
    }
    if generators is not None:
        loaded["generators"] = GeneratorPool(
            [
                HubGenerator(i, backend=backend, models_dir=models_dir)
                for i in generators
            ]
        )
    if similarity is not None:
        loaded["similarity"] = HubSimilarity(
            similarity[0], backend=backend, models_dir=models_dir
        )
    return loaded


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


def _check_repositories(repos: Path | None) -> str:
    if repos is None:
        raise ValueError(
            "the run_sh backend needs the folder with the model repositories: "
            "pass --path-to-model-repositories or set 'path_to_model_repositories' "
            "in the config file"
        )
    if not repos.is_dir():
        raise ValueError(f"the model repositories folder {repos} does not exist")
    return str(repos)


def _check_model_folders(models_dir: str, model_ids: list[str]) -> None:
    missing = [i for i in model_ids if not (Path(models_dir) / i).is_dir()]
    if missing:
        raise ValueError(
            f"no folder for {missing} in {models_dir}; the run_sh backend needs one "
            "folder per model there"
        )


def _to_spec(row: AnnotatorRow, backend: Backend, models_dir: str | None):
    return AnnotatorSpec(
        annotator_id=row.model,
        annotator=build_annotator(
            row.source, row.model, row.column, backend, models_dir
        ),
        cutoff=row.cutoff,
        direction=row.direction,
        weight=1.0 if row.weight is None else row.weight,
        priority=row.priority,
    )
