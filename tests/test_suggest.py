import sys
import types

import pandas as pd
import pytest
from click.testing import CliRunner

from chemsampler.cli import cli
from chemsampler.cli.commands import suggest as suggest_command
from chemsampler.hub import catalog as hub_catalog_module
from chemsampler.hub import precalc
from chemsampler.hub.catalog import CatalogEntry, write_catalog
from chemsampler.hub.precalc import (
    PERCENTILES,
    Pool,
    PrecalcError,
    open_isaura,
    pool_stats,
    require_annotator,
    top_compounds,
)

# --- the query, on real Parquet files ----------------------------------------------------


@pytest.fixture
def con():
    duckdb = pytest.importorskip("duckdb")
    connection = duckdb.connect(":memory:")
    yield connection
    connection.close()


def _write(con, path, values, key="input"):
    """A Parquet file shaped like isaura's: `key`, the SMILES, then the outputs."""
    con.execute(
        f"COPY (SELECT * FROM (VALUES {values}) t(key, {key}, score, count, label)) "
        f"TO '{path}' (FORMAT PARQUET)"
    )
    return str(path)


@pytest.fixture
def chunks(con, tmp_path):
    """Two chunks: a NULL and a NaN score, and `CC` stored in both."""
    first = _write(
        con,
        tmp_path / "chunk_1.parquet",
        "('k1','C',1.0::FLOAT,1,'x'), ('k2','CC',5.0::FLOAT,2,'x'),"
        "('k3','CCC',NULL::FLOAT,NULL::INTEGER,'x'),"
        "('k4','CCCC','NaN'::FLOAT,4,'x')",
    )
    second = _write(
        con,
        tmp_path / "chunk_2.parquet",
        "('k2','CC',5.0::FLOAT,2,'x'), ('k5','CCCCC',9.0::FLOAT,5,'x'),"
        "('k6','O',5.0::FLOAT,6,'x')",
    )
    return [first, second]


def test_the_highest_values_come_first_and_ties_break_on_the_smiles(con, chunks):
    result = top_compounds(con, chunks, "score", "higher", 3)

    assert list(result.columns) == ["smiles", "score"]
    assert result["smiles"].tolist() == ["CCCCC", "CC", "O"]
    assert result["score"].tolist() == [9.0, 5.0, 5.0]


def test_the_lowest_values_come_first_for_lower(con, chunks):
    result = top_compounds(con, chunks, "score", "lower", 2)

    assert result["smiles"].tolist() == ["C", "CC"]
    assert result["score"].tolist() == [1.0, 5.0]


def test_compounds_without_a_value_are_left_out(con, chunks):
    # NaN would otherwise sort above every number, so a failed prediction would
    # be suggested first
    for direction in ("higher", "lower"):
        smiles = top_compounds(con, chunks, "score", direction, 10)["smiles"].tolist()
        assert "CCC" not in smiles and "CCCC" not in smiles
        assert len(smiles) == 4


def test_a_compound_stored_twice_is_listed_once(con, chunks):
    smiles = top_compounds(con, chunks, "score", "higher", 10)["smiles"].tolist()

    assert smiles.count("CC") == 1


def test_a_compound_stored_twice_with_different_values_keeps_the_better_one(
    con, tmp_path
):
    first = _write(con, tmp_path / "a.parquet", "('k1','C',3.0::FLOAT,1,'x')")
    second = _write(con, tmp_path / "b.parquet", "('k1','C',7.0::FLOAT,1,'x')")

    higher = top_compounds(con, [first, second], "score", "higher", 5)
    lower = top_compounds(con, [first, second], "score", "lower", 5)

    assert higher["score"].tolist() == [7.0]
    assert lower["score"].tolist() == [3.0]


def test_duckdbs_progress_bar_is_off_because_it_would_land_in_the_csv(con, chunks):
    con.execute("SET enable_progress_bar = true")

    top_compounds(con, chunks, "score", "higher", 3)

    assert con.execute("SELECT current_setting('enable_progress_bar')").fetchone() == (
        False,
    )


def test_asking_for_more_than_there_is_returns_everything(con, chunks):
    assert len(top_compounds(con, chunks, "score", "higher", 1000)) == 4


def test_an_integer_column_works_too(con, chunks):
    result = top_compounds(con, chunks, "count", "higher", 2)

    assert result["smiles"].tolist() == ["O", "CCCCC"]
    assert result["count"].tolist() == [6, 5]


def test_a_missing_integer_is_left_out_too(con, chunks):
    # there is no isnan for integers, so only the NULL check can drop this one
    smiles = top_compounds(con, chunks, "count", "lower", 1000)["smiles"].tolist()

    assert "CCC" not in smiles
    assert len(smiles) == 5


def test_the_smiles_column_may_be_called_smiles(con, tmp_path):
    path = _write(
        con, tmp_path / "chunk_1.parquet", "('k1','C',1.0::FLOAT,1,'x')", key="smiles"
    )

    result = top_compounds(con, [path], "score", "higher", 1)

    assert list(result.columns) == ["smiles", "score"]


def test_an_unknown_column_says_which_ones_exist(con, chunks):
    with pytest.raises(PrecalcError) as error:
        top_compounds(con, chunks, "nope", "higher", 3)

    assert "no column 'nope'" in str(error.value)
    assert "score, count, label" in str(error.value)  # outputs only, not key/input


def test_a_model_with_thousands_of_columns_gets_a_short_error(con, tmp_path):
    # embedding models have up to 5000 output columns
    outputs = ", ".join(f"{i}::FLOAT AS dim_{i:04d}" for i in range(2048))
    path = tmp_path / "chunk_1.parquet"
    con.execute(
        f"COPY (SELECT 'k1' AS key, 'C' AS input, {outputs}) "
        f"TO '{path}' (FORMAT PARQUET)"
    )

    with pytest.raises(PrecalcError) as error:
        top_compounds(con, [str(path)], "nope", "higher", 3)

    message = str(error.value)
    assert "dim_0000" in message and "dim_0009" in message
    assert "dim_0010" not in message
    assert "2048 columns in all" in message
    assert len(message) < 300


@pytest.mark.parametrize("column", ["input", "key"])
def test_the_columns_isaura_adds_are_not_outputs(con, chunks, column):
    with pytest.raises(PrecalcError, match="no column"):
        top_compounds(con, chunks, column, "higher", 3)


def test_a_text_column_is_refused(con, chunks):
    with pytest.raises(PrecalcError, match="not numbers"):
        top_compounds(con, chunks, "label", "higher", 3)


def test_a_column_name_is_never_run_as_sql(con, chunks):
    with pytest.raises(PrecalcError, match="no column"):
        top_compounds(con, chunks, 'score"; DROP TABLE t; --', "higher", 3)


def test_unreadable_files_are_explained(con, tmp_path):
    with pytest.raises(PrecalcError, match="could not read"):
        top_compounds(con, [str(tmp_path / "missing.parquet")], "score", "higher", 3)


def test_no_files_is_explained(con):
    with pytest.raises(PrecalcError, match="no files"):
        top_compounds(con, [], "score", "higher", 3)


def test_a_direction_or_count_that_makes_no_sense_is_a_bug_not_a_result(con, chunks):
    with pytest.raises(ValueError, match="direction"):
        top_compounds(con, chunks, "score", "up", 3)
    with pytest.raises(ValueError, match="top"):
        top_compounds(con, chunks, "score", "higher", 0)


# --- finding the files, with a stand-in for isaura -----------------------------------------

KEYS = [
    "eos4zfy/v1/tranches/access.json",
    "eos4zfy/v1/tranches/bloom.pkl",
    "eos4zfy/v1/tranches/data/chunk_2.parquet",
    "eos4zfy/v1/tranches/data/chunk_1.parquet",
    "eos4zfy/v2/tranches/data/chunk_1.parquet",
]


@pytest.fixture
def isaura(monkeypatch):
    """A fake `isaura` module: real isaura writes ~/.isaura when imported."""
    pytest.importorskip("botocore")
    state = types.SimpleNamespace(
        keys=list(KEYS), down=False, built=[], list_error=None, chatter=""
    )

    class MinioStore:
        def __init__(self, **kwargs):
            print(state.chatter)  # isaura reports its own problems on stdout
            state.built.append(kwargs)
            if state.down:
                sys.exit(1)

        def list_keys(self, bucket, prefix):
            if state.list_error:
                raise state.list_error
            for key in state.keys:
                if key.startswith(prefix):
                    yield {"Key": key}

    class DuckDBMinio:
        def __init__(self, **kwargs):
            self.con = "connection"

    const = types.SimpleNamespace(
        MINIO_ENDPOINT_CLOUD="http://cloud:9", MINIO_CLOUD_AK="ak", MINIO_CLOUD_SK="sk"
    )
    base = types.ModuleType("isaura.base")
    base.MinioStore, base.DuckDBMinio = MinioStore, DuckDBMinio
    package = types.ModuleType("isaura")
    package.const, package.base = const, base
    monkeypatch.setitem(sys.modules, "isaura", package)
    monkeypatch.setitem(sys.modules, "isaura.base", base)
    state.const = const
    return state


def test_the_files_of_one_version_are_found_in_order(isaura):
    connection, files = open_isaura("eos4zfy", "v1")

    assert connection == "connection"
    assert files == [
        "s3://isaura-public/eos4zfy/v1/tranches/data/chunk_1.parquet",
        "s3://isaura-public/eos4zfy/v1/tranches/data/chunk_2.parquet",
    ]
    assert isaura.built == [
        {"endpoint": "http://cloud:9", "access": "ak", "secret": "sk"}
    ]


def test_a_missing_version_says_which_exist(isaura):
    with pytest.raises(PrecalcError) as error:
        open_isaura("eos4zfy", "v9")

    assert "no version 'v9'" in str(error.value)
    assert "v1, v2" in str(error.value)


def test_a_model_isaura_does_not_hold_is_explained(isaura):
    with pytest.raises(PrecalcError, match="no precalculations for eos0000"):
        open_isaura("eos0000", "v1")


def test_an_unreachable_server_is_explained_not_a_silent_exit(isaura):
    isaura.down = True

    with pytest.raises(PrecalcError, match="cannot reach isaura at http://cloud:9"):
        open_isaura("eos4zfy", "v1")


def test_a_refused_listing_is_explained(isaura):
    from botocore.exceptions import ClientError

    isaura.list_error = ClientError(
        {"Error": {"Code": "403", "Message": "denied"}}, "ListObjectsV2"
    )

    with pytest.raises(PrecalcError, match="could not list isaura-public/eos4zfy"):
        open_isaura("eos4zfy", "v1")


def test_missing_cloud_keys_say_how_to_add_them(isaura):
    isaura.const.MINIO_CLOUD_SK = None

    with pytest.raises(PrecalcError, match="isaura configure --remote"):
        open_isaura("eos4zfy", "v1")
    assert isaura.built == []  # no connection was attempted


def test_isaura_missing_says_how_to_install_it(monkeypatch):
    monkeypatch.setitem(sys.modules, "isaura", None)  # `import isaura` fails

    with pytest.raises(PrecalcError, match=r"pip install chemsampler\[isaura\]"):
        open_isaura("eos4zfy", "v1")


@pytest.mark.parametrize(
    ("model", "version"),
    [("", "v1"), ("eos4zfy/", "v1"), ("eos4zfy", ""), ("eos4zfy", "v1/../v2")],
)
def test_an_odd_model_or_version_never_reaches_the_store(isaura, model, version):
    # an empty model id would list the whole bucket
    with pytest.raises(PrecalcError):
        open_isaura(model, version)

    assert isaura.built == []


def test_what_isaura_prints_goes_to_stderr_not_into_the_csv(isaura, capsys):
    isaura.chatter = "isaura says hello"

    open_isaura("eos4zfy", "v1")

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "isaura says hello" in captured.err


# --- the command ---------------------------------------------------------------------------


def _invoke(*args):
    return CliRunner().invoke(cli, ["suggest", *args])


@pytest.fixture
def hub_catalog(tmp_path, monkeypatch):
    """A small catalog file, so the tests do not depend on the shipped one."""
    path = tmp_path / "hub-catalog.csv"
    write_catalog(
        [
            CatalogEntry("annotators", "eos4zfy", "MAIP", "Annotation", ""),
            CatalogEntry(
                "generators", "eos9p57", "CReM Growth", "Sampling", "Generation"
            ),
            CatalogEntry(
                "similarity",
                "eos8nyg",
                "Reference library",
                "Sampling",
                "Similarity search",
            ),
        ],
        path,
    )
    monkeypatch.setattr(hub_catalog_module, "CATALOG_PATH", path)
    return path


# two compounds the store returned: brackets are what Rich would read as markup
SMILES_A = "CCCNCc1cncc(-c2cnc3[nH]nc(-c4nc5cc(N6CCN(C)CC6)ccc5[nH]4)c3c2)c1"
SMILES_B = "CC1(C)[N+]([O-])=[N+]([O-])C1(C)Br"


@pytest.fixture
def seen(monkeypatch, hub_catalog):
    """Replace the two precalc steps; record what the command asked for."""
    calls = {}

    def fake_open(model, version):
        calls["open"] = (model, version)
        return "connection", ["file"]

    def fake_top(connection, files, column, direction, top):
        calls["top"] = (column, direction, top)
        return pd.DataFrame({"smiles": [SMILES_A, SMILES_B], column: [2.5, 1.0]})

    def fake_pool(connection, files, column):
        calls["pool"] = column
        return Pool(1_355_094, 0.0, 5.0, ((1, 0.1), (50, 2.5), (99.9, 4.99)))

    monkeypatch.setattr(precalc, "open_isaura", fake_open)
    monkeypatch.setattr(precalc, "top_compounds", fake_top)
    monkeypatch.setattr(precalc, "pool_stats", fake_pool)
    return calls


REQUIRED = ["--model", "eos4zfy", "--version", "v2", "--column", "score"]


def test_the_command_prints_only_the_csv_and_asks_for_ten_by_default(seen):
    result = _invoke(*REQUIRED, "--direction", "higher")

    assert result.exit_code == 0, result.output
    assert result.stdout == f"smiles,score\n{SMILES_A},2.5\n{SMILES_B},1.0\n"
    assert seen["open"] == ("eos4zfy", "v2")
    assert seen["top"] == ("score", "higher", 10)
    assert "pool" not in seen  # the extra pass is for a person, not for the CSV


def test_top_can_be_changed(seen):
    _invoke(*REQUIRED, "--direction", "lower", "--top", "3")

    assert seen["top"] == ("score", "lower", 3)


@pytest.mark.parametrize("missing", ["--model", "--version", "--column", "--direction"])
def test_every_option_but_top_is_required(seen, missing):
    args = {**dict(zip(REQUIRED[::2], REQUIRED[1::2])), "--direction": "higher"}
    del args[missing]

    result = _invoke(*(part for pair in args.items() for part in pair))

    assert result.exit_code == 2
    assert missing in result.output
    assert seen == {}


@pytest.mark.parametrize(
    "bad", [["--direction", "up"], ["--direction", "higher", "--top", "0"]]
)
def test_a_bad_direction_or_count_is_a_usage_error(seen, bad):
    result = _invoke(*REQUIRED, *bad)

    assert result.exit_code == 2
    assert seen == {}


def test_a_failure_is_one_red_line_on_stderr_and_nothing_on_stdout(monkeypatch):
    def fail(model, version):
        raise PrecalcError("cannot reach isaura at http://cloud:9")

    monkeypatch.setattr(precalc, "open_isaura", fail)

    result = _invoke(*REQUIRED, "--direction", "higher")

    assert result.exit_code == 1
    assert result.stdout == ""
    assert "cannot reach isaura" in result.stderr


def test_suggest_is_listed_in_the_help_as_a_developer_tool():
    result = CliRunner().invoke(cli, ["--help"], terminal_width=80)

    assert "suggest" in result.output
    assert "developers" in result.output


# --- only annotation models: the ones that give scores -----------------------------------


def test_an_annotation_model_is_accepted(hub_catalog):
    assert require_annotator("eos4zfy").title == "MAIP"


@pytest.mark.parametrize(
    ("model", "kind"),
    [
        ("eos9p57", "Sampling model (Generation)"),
        ("eos8nyg", "Sampling model (Similarity search)"),
    ],
)
def test_a_sampling_model_is_refused_and_says_what_it_is(hub_catalog, model, kind):
    with pytest.raises(PrecalcError) as error:
        require_annotator(model)

    assert kind in str(error.value)
    assert "only works on annotation models" in str(error.value)


def test_a_model_the_catalog_does_not_list_is_refused(hub_catalog):
    # representation (featurization) models are not in the catalog at all
    with pytest.raises(PrecalcError) as error:
        require_annotator("eos4x30")

    assert "not an annotation model" in str(error.value)
    assert "never accepted" in str(error.value)
    assert "chemsampler catalog --update" in str(error.value)  # for a new model


def test_a_malformed_model_id_is_refused_before_the_catalog(hub_catalog):
    with pytest.raises(PrecalcError, match="not a model id"):
        require_annotator("")


def test_an_unreadable_catalog_is_explained(hub_catalog):
    hub_catalog.write_text("garbage")

    with pytest.raises(PrecalcError, match="cannot check eos4zfy against the catalog"):
        require_annotator("eos4zfy")


@pytest.mark.parametrize("model", ["eos9p57", "eos8nyg", "eos4x30"])
def test_the_command_refuses_other_models_before_touching_isaura(seen, model):
    result = _invoke(
        "--model",
        model,
        "--version",
        "v1",
        "--column",
        "score",
        "--direction",
        "higher",
    )

    assert result.exit_code == 1
    assert result.stdout == ""
    assert seen == {}  # isaura was never asked


# --- the table a person sees, and the pool it quotes ---------------------------------------


def test_the_pool_counts_compounds_with_a_value_and_gives_range_and_percentiles(
    con, chunks
):
    # NULL, NaN and the compound stored twice do not count
    quartiles = (25, 50, 75)

    assert pool_stats(con, chunks, "score", quartiles) == Pool(
        4, 1.0, 9.0, ((25, 4.0), (50, 5.0), (75, 6.0))
    )
    assert pool_stats(con, chunks, "count", quartiles) == Pool(
        5, 1.0, 6.0, ((25, 2.0), (50, 4.0), (75, 5.0))
    )


def test_the_default_percentiles_reach_into_the_upper_tail(con, chunks):
    pool = pool_stats(con, chunks, "score")

    levels = [level for level, _ in pool.percentiles]
    values = [value for _, value in pool.percentiles]
    assert levels == list(PERCENTILES) and 99.9 in levels
    assert values == sorted(values)
    assert pool.low <= values[0] and values[-1] <= pool.high


@pytest.mark.parametrize("bad", [(), (0,), (100,), (50, 101), (-5,)])
def test_percentiles_must_be_given_and_inside_zero_to_one_hundred(con, chunks, bad):
    with pytest.raises(ValueError, match="percentiles"):
        pool_stats(con, chunks, "score", bad)


def test_a_pool_with_no_values_is_empty_not_an_error(con, tmp_path):
    path = _write(con, tmp_path / "chunk_1.parquet", "('k1','C',NULL::FLOAT,1,'x')")

    assert pool_stats(con, [path], "score") == Pool(0, 0.0, 0.0)


def test_the_pool_checks_the_column_like_the_query_does(con, chunks):
    with pytest.raises(PrecalcError, match="no column 'nope'"):
        pool_stats(con, chunks, "nope")
    with pytest.raises(PrecalcError, match="not numbers"):
        pool_stats(con, chunks, "label")


@pytest.fixture
def terminal(monkeypatch):
    """Pretend a person is at a terminal, wide enough that nothing wraps."""
    monkeypatch.setattr(suggest_command, "_is_terminal", lambda: True)
    monkeypatch.setenv("COLUMNS", "200")


def test_at_a_terminal_the_result_is_a_table_with_context(seen, terminal):
    result = _invoke(*REQUIRED, "--direction", "higher")

    out = result.stdout
    assert result.exit_code == 0, result.output
    assert "MAIP" in out  # the model's name, from the catalog
    assert "eos4zfy · v2 · score" in out
    assert "2 highest of 1,355,094 scored compounds" in out
    assert "in the Ersilia Reference Library" in out
    assert "> seeds.csv" in out
    assert "smiles,score" not in out  # a table, not the CSV
    assert seen["pool"] == "score"


def test_the_smiles_in_the_table_keep_their_brackets(seen, terminal):
    # Rich would take [nH] for markup and drop it from a bare string
    out = _invoke(*REQUIRED, "--direction", "higher").stdout

    assert SMILES_A in out and SMILES_B in out


def test_the_library_percentiles_sit_between_min_and_max(seen, terminal):
    lines = _invoke(*REQUIRED, "--direction", "higher").stdout.splitlines()

    header = next(i for i, line in enumerate(lines) if "min" in line and "max" in line)
    assert lines[header].split() == ["min", "p1", "p50", "p99.9", "max"]
    assert lines[header + 2].split() == ["0", "0.1", "2.5", "4.99", "5"]


def test_there_is_no_range_column_or_bar_any_more(seen, terminal):
    out = _invoke(*REQUIRED, "--direction", "higher").stdout

    assert "range" not in out
    assert "█" not in out and "░" not in out


def test_lower_says_lowest(seen, terminal):
    out = _invoke(*REQUIRED, "--direction", "lower").stdout

    assert "2 lowest of" in out


def test_a_column_without_values_says_so_and_skips_the_extra_pass(
    seen, terminal, monkeypatch
):
    monkeypatch.setattr(
        precalc,
        "top_compounds",
        lambda *args: pd.DataFrame({"smiles": [], "score": []}),
    )

    result = _invoke(*REQUIRED, "--direction", "higher")

    assert result.exit_code == 0, result.output
    assert "No compound has a value for this column" in result.stdout
    assert "pool" not in seen
