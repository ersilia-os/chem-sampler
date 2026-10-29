import pytest

from chemsampler.config.schema import check_header, parse_annotator, parse_cutoff
from chemsampler.models.spec import AnnotatorSpec


def _row(**overrides):
    row = {
        "source": "ersilia",
        "model": "eos4zfy",
        "cutoff": "0.5",
        "direction": "higher",
        "priority": "primary",
    }
    row.update(overrides)
    return row


# --- cutoff and direction ---------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("0.3", 0.3),
        (0.3, 0.3),
        (5, 5.0),
        ("1e-3", 0.001),
        ("-2", -2.0),
        ("-inf", float("-inf")),
    ],
)
def test_single_cutoffs(raw, expected):
    assert parse_cutoff("mw", "higher", raw) == expected


@pytest.mark.parametrize(
    "raw, expected",
    [("250-450", (250.0, 450.0)), (" 0.2 - 0.8 ", (0.2, 0.8)), ("-5--1", (-5.0, -1.0))],
)
def test_range_cutoffs(raw, expected):
    assert parse_cutoff("mw", "range", raw) == expected


def test_range_direction_needs_a_range_and_says_how_to_write_one():
    with pytest.raises(ValueError) as error:
        parse_cutoff("mw", "range", "250")

    message = str(error.value)
    assert "direction is 'range'" in message
    assert "'min-max'" in message and "250-450" in message
    assert "'higher' or 'lower'" in message


def test_one_sided_direction_with_a_range_says_to_use_range():
    with pytest.raises(ValueError, match="is a range, so direction must be 'range'"):
        parse_cutoff("mw", "higher", "250-450")


def test_a_range_needs_min_below_max():
    with pytest.raises(ValueError, match="min must be lower than max"):
        parse_cutoff("mw", "range", "450-250")


@pytest.mark.parametrize("raw", ["", "  ", None, True])
def test_a_cutoff_is_required(raw):
    with pytest.raises(ValueError, match="cutoff is required"):
        parse_cutoff("mw", "higher", raw)


def test_an_unreadable_cutoff_is_rejected():
    with pytest.raises(ValueError, match="invalid cutoff"):
        parse_cutoff("mw", "higher", "high")


# --- a whole annotator entry ------------------------------------------------------


def test_parses_a_primary_ersilia_annotator():
    row = parse_annotator(_row(weight="2", column="score"), "annotators.csv row 1")

    assert (row.source, row.model, row.cutoff, row.direction) == (
        "ersilia",
        "eos4zfy",
        0.5,
        "higher",
    )
    assert (row.priority, row.weight, row.column) == ("primary", 2.0, "score")


def test_parses_a_secondary_rdkit_guard_rail():
    row = parse_annotator(
        _row(
            source="rdkit",
            model="mw",
            cutoff="250-450",
            direction="range",
            priority="secondary",
        ),
        "annotators.csv row 2",
    )

    assert row.cutoff == (250.0, 450.0)
    assert row.weight is None and row.column is None


def test_empty_cells_and_whitespace_count_as_missing_or_are_stripped():
    row = parse_annotator(_row(model=" eos4zfy ", weight="", column=" "), "row")

    assert (row.model, row.weight, row.column) == ("eos4zfy", None, None)


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"source": "pubchem"}, "source must be one of"),
        ({"direction": "up"}, "direction must be one of"),
        ({"priority": ""}, "missing"),
        ({"priority": "main"}, "priority must be one of"),
        ({"weight": "-1"}, "weight must be >= 0"),
        ({"weight": "heavy"}, "invalid weight"),
        ({"cutoff": "0.5-0.9", "direction": "range"}, "only allowed for secondary"),
        ({"priority": "secondary", "weight": "2"}, "weight only applies to primary"),
        ({"source": "rdkit", "model": "logp"}, "unknown rdkit model"),
        ({"source": "rdkit", "model": "qed", "column": "x"}, "column only applies"),
        ({"colour": "red"}, "unknown field"),
    ],
)
def test_invalid_entries_are_explained(overrides, message):
    with pytest.raises(ValueError, match=message):
        parse_annotator(_row(**overrides), "annotators.csv row 3")


def test_errors_say_where_they_are():
    with pytest.raises(ValueError, match=r"annotators\.csv row 3 \(eos4zfy\)"):
        parse_annotator(_row(priority="main"), "annotators.csv row 3")


# --- the spec keeps the same invariants for Python callers -------------------------


def test_spec_rejects_a_range_on_a_primary():
    with pytest.raises(ValueError, match="only allowed for secondary"):
        AnnotatorSpec("mw", object(), (250.0, 450.0), "range")


def test_spec_rejects_a_range_without_range_direction_and_vice_versa():
    with pytest.raises(ValueError, match="direction must be 'range'"):
        AnnotatorSpec("mw", object(), (250.0, 450.0), "higher", priority="secondary")
    with pytest.raises(ValueError, match="needs a \\(min, max\\) cutoff"):
        AnnotatorSpec("mw", object(), 250.0, "range", priority="secondary")


def test_spec_rejects_a_bad_range_and_an_unknown_priority():
    with pytest.raises(ValueError, match="min < max"):
        AnnotatorSpec("mw", object(), (450.0, 250.0), "range", priority="secondary")
    with pytest.raises(ValueError, match="priority must be"):
        AnnotatorSpec("mw", object(), 1.0, "higher", priority="main")


# --- CSV headers ------------------------------------------------------------------


def test_header_must_match():
    check_header("g.csv", ["generators"], ["generators"])
    check_header(
        "a.csv",
        ["source", "model", "cutoff", "direction", "priority", "weight"],
        ["source", "model", "cutoff", "direction", "priority"],
        ["weight", "column"],
    )


@pytest.mark.parametrize(
    "header, message",
    [
        (["generator_id"], "missing"),
        (["generators", "extra"], "unknown"),
        (["generators", "generators"], "repeated"),
        (None, "missing"),
    ],
)
def test_bad_headers_are_explained(header, message):
    with pytest.raises(ValueError, match=message):
        check_header("generators.csv", header, ["generators"])
