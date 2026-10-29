"""Top-N compounds from isaura's precalculations, for `chemsampler suggest`.

isaura keeps the outputs of Hub models for the Ersilia Reference Library as
Parquet files on an S3-compatible store. Everything that touches isaura is in
`open_isaura`; `top_compounds` needs only a DuckDB connection and file URIs, so it
runs on local files too. Nothing here is shared with the other commands.
"""

import contextlib
import re
import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .catalog import CatalogEntry

BUCKET = "isaura-public"

#: Columns isaura adds itself; the rest are the model's outputs.
KEY_COLUMNS = ("key", "input", "smiles")
SMILES_COLUMNS = ("input", "smiles")
DIRECTIONS = ("higher", "lower")

_MODEL_ID = re.compile(r"eos\d[a-z0-9]{3}")
_VERSION = re.compile(r"[A-Za-z0-9._-]+")
_NUMBERS = {
    "TINYINT",
    "SMALLINT",
    "INTEGER",
    "BIGINT",
    "HUGEINT",
    "UTINYINT",
    "USMALLINT",
    "UINTEGER",
    "UBIGINT",
    "UHUGEINT",
    "FLOAT",
    "DOUBLE",
    "DECIMAL",
}
_FLOATS = {"FLOAT", "DOUBLE"}


class PrecalcError(Exception):
    """Something the user can act on: it is printed as one line, without a traceback."""


#: Where along the library's values `pool_stats` reports by default. The upper tail is
#: dense on purpose: suggestions come from the very top, far beyond the 99th.
PERCENTILES = (1, 5, 25, 50, 75, 95, 99, 99.9)


@dataclass(frozen=True)
class Pool:
    """
    The compounds that have a value for a column.

    `percentiles` pairs each requested percentile with the value at that point;
    it is empty when no compound has a value.
    """

    count: int
    low: float
    high: float
    percentiles: tuple[tuple[float, float], ...] = ()


def _quote(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _source(files: list[str]) -> str:
    uris = ", ".join("'" + uri.replace("'", "''") + "'" for uri in files)
    return f"read_parquet([{uris}])"


def _execute(con, sql: str):
    import duckdb  # a DuckDB connection is in hand, so the extra is installed

    try:
        return con.execute(sql)
    except duckdb.Error as e:
        raise PrecalcError(f"could not read the precalculations: {e}") from e


def _listed(names: list[str], limit: int = 10) -> str:
    """The names, or the first few and a count: some models have thousands of columns."""
    if len(names) <= limit:
        return ", ".join(names) or "none"
    return f"{', '.join(names[:limit])}, ... ({len(names)} columns in all)"


def _base_type(duckdb_type: str) -> str:
    """`DECIMAL(18,3)` -> `DECIMAL`."""
    return duckdb_type.split("(")[0].strip().upper()


def _prepare(con, files: list[str], column: str) -> tuple[str, str, str, str]:
    """
    Check `column` against the files and build the SQL pieces both queries share.

    Returns `(source, key, value, present)`: the `read_parquet` expression, the quoted
    SMILES column, the quoted value column, and the condition that keeps only
    compounds that have a value.
    """
    if not files:
        raise PrecalcError("there are no files to read")

    # DuckDB draws a progress bar on stdout during a long query (a version split into
    # many chunks gets there), and that would end up inside the CSV
    _execute(con, "SET enable_progress_bar = false")
    source = _source(files)
    rows = _execute(con, f"DESCRIBE SELECT * FROM {source}").fetchall()
    types = {name: kind for name, kind, *_ in rows}
    key = next((name for name in SMILES_COLUMNS if name in types), None)
    if key is None:
        raise PrecalcError("the files have no `input` or `smiles` column")
    outputs = [name for name in types if name not in KEY_COLUMNS]
    # the name goes into the SQL, so it must be one the files really have
    if column not in outputs:
        raise PrecalcError(f"no column {column!r}; this model has: {_listed(outputs)}")
    kind = _base_type(types[column])
    if kind not in _NUMBERS:
        raise PrecalcError(f"column {column!r} holds {types[column]}, not numbers")

    value = _quote(column)
    present = f"{value} IS NOT NULL"
    if kind in _FLOATS:
        present += f" AND NOT isnan({value})"
    return source, _quote(key), value, present


def top_compounds(
    con, files: list[str], column: str, direction: str, top: int
) -> pd.DataFrame:
    """
    The compounds with the highest or lowest value of one output column.

    Parameters
    ----------
    con
        A DuckDB connection that can read `files`.
    files : list of str
        Parquet files (paths or S3 URIs) holding one model version.
    column : str
        An output column of the model.
    direction : {"higher", "lower"}
        Which end of the column to take.
    top : int
        How many compounds to return.

    Returns
    -------
    pandas.DataFrame
        Columns `smiles` and `column`, best first, one row per compound. Compounds
        without a value (NULL or NaN) are left out: NaN would sort above every number.

    Raises
    ------
    PrecalcError
        If the files cannot be read, or `column` is not a numeric output column.
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, not {direction!r}")
    if top < 1:
        raise ValueError(f"top must be at least 1, not {top}")
    source, key, value, present = _prepare(con, files, column)
    pick, order = ("max", "DESC") if direction == "higher" else ("min", "ASC")
    # one row per compound: some versions store a few compounds twice
    sql = (
        f"SELECT smiles, score AS {value} FROM ("
        f"SELECT {key} AS smiles, {pick}({value}) AS score "
        f"FROM {source} WHERE {present} GROUP BY {key}) "
        f"ORDER BY score {order}, smiles LIMIT {int(top)}"
    )
    return _execute(con, sql).df()


def pool_stats(
    con, files: list[str], column: str, percentiles: tuple[float, ...] = PERCENTILES
) -> Pool:
    """
    How many compounds have a value for `column`, their range, and some percentiles.

    Each compound counts once. It is one more pass over the column, so it is for
    showing context next to `top_compounds`, not for the CSV.

    Parameters
    ----------
    con
        A DuckDB connection that can read `files`.
    files : list of str
        Parquet files holding one model version.
    column : str
        An output column of the model.
    percentiles : tuple of float
        Percentiles to report, each strictly between 0 and 100.

    Raises
    ------
    PrecalcError
        If the files cannot be read, or `column` is not a numeric output column.
    """
    if not percentiles or not all(0 < p < 100 for p in percentiles):
        raise ValueError(
            f"percentiles must be given, between 0 and 100 exclusive: {percentiles}"
        )
    source, key, value, present = _prepare(con, files, column)
    # a compound stored twice has the same value both times, so any pick will do
    levels = ", ".join(str(round(p / 100, 6)) for p in percentiles)
    sql = (
        f"SELECT count(*), min(score), max(score), quantile_cont(score, [{levels}]) "
        f"FROM (SELECT max({value}) AS score FROM {source} "
        f"WHERE {present} GROUP BY {key})"
    )
    count, low, high, cuts = _execute(con, sql).fetchone()
    if not count:
        return Pool(0, 0.0, 0.0)
    return Pool(
        int(count),
        float(low),
        float(high),
        tuple((p, float(cut)) for p, cut in zip(percentiles, cuts)),
    )


def _check_model(model: str) -> None:
    # an empty model would list the whole bucket, so be strict
    if not _MODEL_ID.fullmatch(model):
        raise PrecalcError(f"{model!r} is not a model id like eos4zfy")


def _check_names(model: str, version: str) -> None:
    _check_model(model)
    if not _VERSION.fullmatch(version):
        raise PrecalcError(f"{version!r} is not a version like v1")


def require_annotator(model: str) -> "CatalogEntry":
    """
    Accept only annotation models, the ones whose outputs are scores.

    A top-N only means something for a score. The columns of a representation
    (featurization) model are embedding dimensions, and a sampling model
    (generator, similarity search) has no score at all. The bundled catalog gives
    the Hub's label for each model, so this needs no network and runs before isaura
    is touched. A model the catalog does not list is refused too.

    Returns
    -------
    CatalogEntry
        The model's catalog entry, whose title is the model's name.

    Raises
    ------
    PrecalcError
        If `model` is not a catalogued annotation model.
    """
    from .catalog import read_catalog

    _check_model(model)
    try:
        entry = read_catalog().by_id().get(model)
    except ValueError as e:
        raise PrecalcError(f"cannot check {model} against the catalog: {e}") from e
    if entry is None:
        raise PrecalcError(
            f"{model} is not an annotation model in the bundled catalog, so it has no "
            "scores to rank (representation and sampling models are never accepted). "
            "If it is a new annotation model, run `chemsampler catalog --update`"
        )
    if entry.role != "annotators":
        raise PrecalcError(
            f"{model} ({entry.title}) is a {entry.task} model"
            f"{f' ({entry.subtask})' if entry.subtask else ''}; "
            "suggest only works on annotation models, which give scores"
        )
    return entry


def open_isaura(model: str, version: str, bucket: str = BUCKET) -> tuple:
    """
    Connect to isaura's cloud store and find the files of one model version.

    The endpoint and keys are the ones `isaura configure --remote` saved. isaura is
    imported here and nowhere else: importing it writes `~/.isaura/` if absent.

    Parameters
    ----------
    model : str
        Ersilia model id, e.g. `eos4zfy`.
    version : str
        The stored version, e.g. `v2`. There is deliberately no default: versions
        of one model can differ.
    bucket : str
        The isaura bucket to read.

    Returns
    -------
    tuple
        `(connection, files)` for `top_compounds`.

    Raises
    ------
    PrecalcError
        If isaura is missing or unreachable, or holds nothing for that version.
    """
    _check_names(model, version)
    try:
        from botocore.exceptions import BotoCoreError, ClientError
        from isaura import const
        from isaura.base import DuckDBMinio, MinioStore
    except ImportError as e:
        raise PrecalcError(
            f"isaura could not be imported ({e}); "
            "install it with `pip install chemsampler[isaura]`"
        ) from e

    endpoint = const.MINIO_ENDPOINT_CLOUD
    access, secret = const.MINIO_CLOUD_AK, const.MINIO_CLOUD_SK
    if not (access and secret):
        raise PrecalcError("no cloud keys for isaura; run `isaura configure --remote`")

    # isaura reports its own problems on stdout, which has to stay a clean CSV
    with contextlib.redirect_stdout(sys.stderr):
        try:
            store = MinioStore(endpoint=endpoint, access=access, secret=secret)
        except SystemExit as e:  # isaura exits when the server does not answer
            raise PrecalcError(f"cannot reach isaura at {endpoint}") from e
        try:
            keys = [obj["Key"] for obj in store.list_keys(bucket, f"{model}/")]
        except (BotoCoreError, ClientError) as e:
            raise PrecalcError(f"could not list {bucket}/{model}: {e}") from e
        connection = DuckDBMinio(endpoint=endpoint, access=access, secret=secret).con

    prefix = f"{model}/{version}/"
    files = sorted(
        f"s3://{bucket}/{key}"
        for key in keys
        if key.startswith(prefix) and key.endswith(".parquet") and "/chunk_" in key
    )
    if not files:
        versions = sorted({key.split("/")[1] for key in keys if key.count("/") > 1})
        if versions:
            raise PrecalcError(
                f"{model} has no version {version!r} in {bucket}; "
                f"it has: {', '.join(versions)}"
            )
        raise PrecalcError(f"{bucket} has no precalculations for {model}")
    return connection, files
