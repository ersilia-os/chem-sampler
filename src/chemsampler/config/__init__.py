import csv
from pathlib import Path

from ..hub.client import Backend
from ..models.annotator import QEDAnnotator
from ..models.chembl import ChemblSampler
from ..models.generator import GeneratorPool, HubGenerator
from ..models.hub_annotator import HubAnnotator
from ..models.spec import AnnotatorSpec

#: Final fallback for `load_generators()`, used only when neither an explicit
#: path nor a ./generators.csv in the cwd is available. Deliberately excludes
#: "chembl": the ChemblSampler baseline must only ever be opted into by a
#: user-supplied CSV.
_DEFAULT_GENERATORS_PATH = Path(__file__).parent / "generators.csv"


def load_generators(
    path: str | None = None, backend: Backend = "run_sh"
) -> GeneratorPool:
    """
    Build a pool of generators from a CSV of generator ids.

    Parameters
    ----------
    path : str, optional
        CSV with a `generator_id` column, one row per generator. Each id is
        either an Ersilia Hub model id, or the literal "chembl" to opt into
        `ChemblSampler`. Resolved in order: `path` if given; else
        `./generators.csv` in the current working directory, if one exists;
        else the generators validated in ersilia-os/ersilia#1919 (never
        includes "chembl" - must be listed explicitly in a user-supplied CSV).
    backend : {"ersilia", "run_sh"}, optional
        Passed through to every `HubGenerator` built from this CSV, by default
        "run_sh". Has no effect on a "chembl" row (`ChemblSampler` has no
        backend concept).

    Returns
    -------
    GeneratorPool
        Pool of the generators listed in the CSV.
    """
    cwd_path = Path("generators.csv")
    if path is not None:
        resolved = path
    elif cwd_path.exists():
        resolved = cwd_path
    else:
        resolved = _DEFAULT_GENERATORS_PATH
    rows = _read_csv(resolved, required_columns=("generator_id",))
    return GeneratorPool(
        [_build_generator(row["generator_id"], backend) for row in rows]
    )


def load_annotators(
    path: str | None = None, backend: Backend = "run_sh"
) -> list[AnnotatorSpec]:
    """
    Build annotator specs from a CSV of annotator ids, cutoffs and directions.

    Parameters
    ----------
    path : str, optional
        CSV with columns `annotator_id, cutoff, direction, column`.
        `annotator_id` is either "qed" or an Ersilia Hub model id. `cutoff` and
        `direction` ("higher" or "lower") are required for every row - every
        annotator is uniform, there is no default. Row order is meaningful: it's
        the priority order `hill_climb` uses in `mode="sequential"`. `column` is
        optional and is passed to `HubAnnotator` for models with more than one
        numeric output. Resolved in order: `path` if given; else
        `./annotators.csv` in the current working directory, if one exists.
        Unlike `load_generators`, there is no further built-in default - see
        Raises.
    backend : {"ersilia", "run_sh"}, optional
        Passed through to every `HubAnnotator` built from this CSV, by default
        "run_sh". Has no effect on a "qed" row (`QEDAnnotator` has no backend
        concept).

    Returns
    -------
    list[AnnotatorSpec]
        One spec per row, in file order.

    Raises
    ------
    FileNotFoundError
        If `path` is None and no `./annotators.csv` exists in the current
        working directory. Annotators have no built-in default the way
        generators do - cutoffs and directions are always user-specific.
    ValueError
        If a required column is missing, a cutoff or direction is invalid, or
        the annotator ids are not unique.
    """
    cwd_path = Path("annotators.csv")
    if path is not None:
        resolved = path
    elif cwd_path.exists():
        resolved = cwd_path
    else:
        raise FileNotFoundError(
            "no annotators path given and no ./annotators.csv in the current "
            "directory; annotators have no built-in default, unlike "
            "generators (cutoffs and directions are always user-specific)"
        )

    rows = _read_csv(resolved, required_columns=("annotator_id", "cutoff", "direction"))

    specs = [
        AnnotatorSpec(
            annotator_id=row["annotator_id"],
            annotator=build_annotator(
                row["annotator_id"], row.get("column") or None, backend
            ),
            cutoff=_parse_cutoff(row["annotator_id"], row["cutoff"]),
            direction=row["direction"],
        )
        for row in rows
    ]

    ids = [spec.annotator_id for spec in specs]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{resolved}: duplicate annotator_id")

    return specs


def _build_generator(generator_id: str, backend: Backend = "run_sh"):
    """Resolve a generator id to a generator instance."""
    if generator_id == "chembl":
        return ChemblSampler()
    return HubGenerator(generator_id, backend=backend)


def build_annotator(
    annotator_id: str, column: str | None = None, backend: Backend = "run_sh"
):
    """
    Resolve an annotator id to an annotator instance.

    Parameters
    ----------
    annotator_id : str
        Either "qed", or an Ersilia Hub model id.
    column : str, optional
        Output column to score on, passed through to `HubAnnotator`. Ignored
        for "qed". By default `None` (the model's single numeric output;
        required if it returns more than one).
    backend : {"ersilia", "run_sh"}, optional
        Passed through to `HubAnnotator`, by default "run_sh". Ignored for "qed".

    Returns
    -------
    QEDAnnotator or HubAnnotator
        `QEDAnnotator` if `annotator_id == "qed"`, otherwise `HubAnnotator`.
    """
    if annotator_id == "qed":
        return QEDAnnotator()
    return HubAnnotator(annotator_id, column=column, backend=backend)


def _parse_cutoff(annotator_id: str, value: str) -> float:
    """A cutoff is required for every annotator - fail loud rather than guess."""
    if not value:
        raise ValueError(f"{annotator_id}: cutoff is required")
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"{annotator_id}: invalid cutoff {value!r}") from None


def _read_csv(path, required_columns: tuple[str, ...]) -> list[dict]:
    """Read a CSV as a list of row dicts, failing loudly on missing columns."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in required_columns if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} is missing required column(s): {missing}")
        return list(reader)
