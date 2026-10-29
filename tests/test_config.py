import csv
import json

import pytest

from chemsampler import config
from chemsampler.config import (
    build_annotator,
    load_annotators,
    load_inputs,
    plan_fetch,
)
from chemsampler.hub.envs import env_path, write_marker
from chemsampler.models.annotator import MolecularWeightAnnotator, QEDAnnotator
from chemsampler.models.hub_annotator import HubAnnotator
from chemsampler.models.similarity import HubSimilarity
from chemsampler.workspace import create_workspace

ANNOTATOR_HEADER = ["source", "model", "cutoff", "direction", "priority", "weight"]
ANNOTATOR_ROWS = [
    ["ersilia", "eos4zfy", "0.5", "higher", "primary", "2"],
    ["rdkit", "mw", "250-450", "range", "secondary", ""],
    ["rdkit", "qed", "0.3", "higher", "secondary", ""],
]


@pytest.fixture(autouse=True)
def hub(monkeypatch):
    """Record what would be sent to the Hub instead of querying it."""
    calls = []
    monkeypatch.setattr(config, "validate_models", lambda **kw: calls.append(kw))
    return calls


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


@pytest.fixture
def inputs_dir(tmp_path):
    folder = tmp_path / "inputs"
    folder.mkdir()
    _write_csv(folder / "similarity.csv", ["similarity"], [["eos8nyg"]])
    _write_csv(folder / "generators.csv", ["generators"], [["eos9taz"], ["eos80sj"]])
    _write_csv(folder / "annotators.csv", ANNOTATOR_HEADER, ANNOTATOR_ROWS)
    return folder


def _load(inputs_dir, **kwargs):
    kwargs.setdefault("backend", "ersilia")
    return load_inputs(path_to_models=str(inputs_dir), **kwargs)


# --- the three CSVs -----------------------------------------------------------------


def test_loads_the_three_inputs_from_csvs(inputs_dir):
    inputs = _load(inputs_dir)

    assert isinstance(inputs.similarity, HubSimilarity)
    assert inputs.similarity.model_id == "eos8nyg"
    assert [g.model_id for g in inputs.generators.generators] == ["eos9taz", "eos80sj"]
    specs = {s.annotator_id: s for s in inputs.annotators}
    assert list(specs) == ["eos4zfy", "mw", "qed"]
    assert isinstance(specs["eos4zfy"].annotator, HubAnnotator)
    assert (specs["eos4zfy"].priority, specs["eos4zfy"].weight) == ("primary", 2.0)
    assert isinstance(specs["mw"].annotator, MolecularWeightAnnotator)
    assert isinstance(specs["qed"].annotator, QEDAnnotator)
    assert (specs["mw"].cutoff, specs["mw"].direction) == ((250.0, 450.0), "range")
    assert specs["mw"].priority == "secondary"


def test_every_hub_model_is_checked_but_in_house_annotators_are_not(inputs_dir, hub):
    _load(inputs_dir)

    assert hub == [
        {
            "generators": ["eos9taz", "eos80sj"],
            "similarity": ["eos8nyg"],
            "annotators": ["eos4zfy"],
        }
    ]


def test_annotators_are_required_and_the_error_lists_every_way(inputs_dir):
    (inputs_dir / "annotators.csv").unlink()

    with pytest.raises(ValueError) as error:
        _load(inputs_dir)

    message = str(error.value)
    assert "annotators is required" in message
    for way in (
        "--annotators FILE",
        "config file",
        "--path-to-models",
        "working directory",
    ):
        assert way in message


def test_similarity_and_generators_are_each_optional(inputs_dir, hub):
    (inputs_dir / "similarity.csv").unlink()
    without_similarity = _load(inputs_dir)
    (inputs_dir / "similarity.csv").write_text("similarity\neos8nyg\n")
    (inputs_dir / "generators.csv").unlink()
    without_generators = _load(inputs_dir)

    assert without_similarity.similarity is None
    assert without_similarity.generators is not None
    assert without_generators.generators is None
    assert without_generators.similarity.model_id == "eos8nyg"
    assert hub[0]["similarity"] == [] and hub[1]["generators"] == []


def test_at_least_one_of_similarity_and_generators_is_needed(inputs_dir):
    (inputs_dir / "similarity.csv").unlink()
    (inputs_dir / "generators.csv").unlink()

    with pytest.raises(ValueError, match="at least one of similarity and generators"):
        _load(inputs_dir)


def test_without_any_source_the_error_lists_the_ways(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError, match="config file.*--path-to-models"):
        load_inputs(backend="ersilia")


# --- flags and the working directory ------------------------------------------------


def test_a_flag_overrides_the_config_and_the_folder(inputs_dir, tmp_path):
    other = tmp_path / "other.csv"
    _write_csv(other, ["generators"], [["eos8vud"]])

    inputs = _load(inputs_dir, generators=str(other))

    assert [g.model_id for g in inputs.generators.generators] == ["eos8vud"]
    assert inputs.similarity.model_id == "eos8nyg"  # the other inputs are unaffected


def test_a_flag_overrides_an_inline_config_value(tmp_path):
    other = tmp_path / "annotators.csv"
    _write_csv(other, ANNOTATOR_HEADER, ANNOTATOR_ROWS[:1])
    path = _write_config(tmp_path / "config.json", **INLINE)

    inputs = load_inputs(config=path, annotators=str(other), backend="ersilia")

    assert [s.annotator_id for s in inputs.annotators] == ["eos4zfy"]


def test_a_flag_file_needs_the_right_header(inputs_dir, tmp_path):
    other = tmp_path / "other.csv"
    _write_csv(other, ["generator_id"], [["eos8vud"]])

    with pytest.raises(ValueError, match="bad header"):
        _load(inputs_dir, generators=str(other))


def test_inputs_are_found_in_the_working_directory(inputs_dir, monkeypatch):
    monkeypatch.chdir(inputs_dir)

    inputs = load_inputs(backend="ersilia")

    assert inputs.similarity.model_id == "eos8nyg"
    assert [s.annotator_id for s in inputs.annotators] == ["eos4zfy", "mw", "qed"]


def test_a_named_folder_never_falls_back_to_the_working_directory(
    inputs_dir, tmp_path, monkeypatch
):
    folder = tmp_path / "folder"
    folder.mkdir()
    (inputs_dir / "similarity.csv").rename(tmp_path / "similarity.csv")
    (inputs_dir / "generators.csv").rename(folder / "generators.csv")
    (inputs_dir / "annotators.csv").rename(folder / "annotators.csv")
    monkeypatch.chdir(tmp_path)  # holds a similarity.csv that must be ignored

    inputs = _load(folder)

    assert inputs.similarity is None


def test_path_to_models_must_be_a_folder(tmp_path):
    with pytest.raises(ValueError, match="is not a folder"):
        load_inputs(path_to_models=str(tmp_path / "nowhere"), backend="ersilia")


def test_csv_headers_must_be_the_input_names(inputs_dir):
    _write_csv(inputs_dir / "generators.csv", ["generator_id"], [["eos9taz"]])

    with pytest.raises(ValueError, match="bad header"):
        _load(inputs_dir)


def test_blank_rows_are_skipped_and_repeats_rejected(inputs_dir):
    _write_csv(
        inputs_dir / "generators.csv", ["generators"], [["eos9taz"], [""], [" "]]
    )
    assert [g.model_id for g in _load(inputs_dir).generators.generators] == ["eos9taz"]

    _write_csv(
        inputs_dir / "generators.csv", ["generators"], [["eos9taz"], ["eos9taz"]]
    )
    with pytest.raises(ValueError, match="more than once"):
        _load(inputs_dir)


def test_an_empty_list_is_rejected(inputs_dir):
    _write_csv(inputs_dir / "generators.csv", ["generators"], [])

    with pytest.raises(ValueError, match="at least one model"):
        _load(inputs_dir)


def test_only_one_similarity_model_for_now(inputs_dir):
    _write_csv(
        inputs_dir / "similarity.csv", ["similarity"], [["eos8nyg"], ["eos2a9n"]]
    )

    with pytest.raises(ValueError, match="only one similarity model"):
        _load(inputs_dir)


def test_annotator_errors_say_the_file_and_row(inputs_dir):
    rows = [ANNOTATOR_ROWS[0], ["rdkit", "mw", "250", "range", "secondary", ""]]
    _write_csv(inputs_dir / "annotators.csv", ANNOTATOR_HEADER, rows)

    with pytest.raises(ValueError, match=r"annotators\.csv row 2 \(mw\).*min-max"):
        _load(inputs_dir)


def test_a_run_needs_a_primary_annotator(inputs_dir):
    _write_csv(inputs_dir / "annotators.csv", ANNOTATOR_HEADER, ANNOTATOR_ROWS[1:])

    with pytest.raises(ValueError, match="at least one primary"):
        _load(inputs_dir)


def test_an_annotator_cannot_be_listed_twice(inputs_dir):
    rows = [ANNOTATOR_ROWS[0], ANNOTATOR_ROWS[0]]
    _write_csv(inputs_dir / "annotators.csv", ANNOTATOR_HEADER, rows)

    with pytest.raises(ValueError, match="more than once"):
        _load(inputs_dir)


def test_weight_defaults_to_one_and_column_reaches_the_hub_annotator(inputs_dir):
    header = ["source", "model", "cutoff", "direction", "priority", "column"]
    rows = [["ersilia", "eos4zfy", "0.5", "higher", "primary", "score"]]
    _write_csv(inputs_dir / "annotators.csv", header, rows)

    spec = _load(inputs_dir).annotators[0]

    assert spec.weight == 1.0
    assert spec.annotator.column == "score"


# --- the config file ----------------------------------------------------------------

INLINE = {
    "similarity": ["eos8nyg"],
    "generators": ["eos9taz"],
    "annotators": [
        {
            "source": "ersilia",
            "model": "eos4zfy",
            "cutoff": 0.5,
            "direction": "higher",
            "priority": "primary",
        },
        {
            "source": "rdkit",
            "model": "mw",
            "cutoff": "250-450",
            "direction": "range",
            "priority": "secondary",
        },
    ],
}


def _write_config(path, **keys):
    path.write_text(json.dumps(keys))
    return str(path)


def test_everything_can_be_given_inline_in_the_config(tmp_path):
    path = _write_config(tmp_path / "config.json", **INLINE)

    inputs = load_inputs(config=path, backend="ersilia")

    assert inputs.similarity.model_id == "eos8nyg"
    assert [s.annotator_id for s in inputs.annotators] == ["eos4zfy", "mw"]
    assert inputs.annotators[1].cutoff == (250.0, 450.0)


def test_config_paths_are_relative_to_the_config_file_not_the_cwd(
    tmp_path, monkeypatch
):
    config_dir = tmp_path / "project"
    config_dir.mkdir()
    (config_dir / "models").mkdir()
    _write_csv(config_dir / "models" / "generators.csv", ["generators"], [["eos9taz"]])
    path = _write_config(
        config_dir / "config.json",
        path_to_models="models",
        similarity=INLINE["similarity"],
        annotators=INLINE["annotators"],
    )
    monkeypatch.chdir(tmp_path)

    inputs = load_inputs(config=path, backend="ersilia")

    assert [g.model_id for g in inputs.generators.generators] == ["eos9taz"]


def test_flags_override_the_config_paths(tmp_path, inputs_dir):
    path = _write_config(tmp_path / "config.json", path_to_models="nowhere")

    inputs = load_inputs(config=path, path_to_models=str(inputs_dir), backend="ersilia")

    assert inputs.similarity.model_id == "eos8nyg"


def test_an_input_in_both_the_config_and_a_csv_is_rejected(tmp_path, inputs_dir):
    path = _write_config(
        tmp_path / "config.json", path_to_models=str(inputs_dir), generators=["eos9taz"]
    )

    with pytest.raises(ValueError, match="generators is given twice"):
        load_inputs(config=path, backend="ersilia")


@pytest.mark.parametrize(
    "keys, message",
    [
        ({"colour": 1}, "unknown key"),
        ({"path_to_models": 3}, "must be a path"),
        ({"generators": "eos9taz"}, "list of model ids"),
        ({"annotators": ["mw"]}, "list of objects"),
    ],
)
def test_malformed_configs_are_explained(tmp_path, keys, message):
    path = _write_config(tmp_path / "config.json", **keys)

    with pytest.raises(ValueError, match=message):
        load_inputs(config=path, backend="ersilia")


def test_config_must_be_a_json_object(tmp_path):
    path = tmp_path / "config.json"
    path.write_text("[1, 2]")
    with pytest.raises(ValueError, match="expected a JSON object"):
        load_inputs(config=str(path), backend="ersilia")

    path.write_text("{oops")
    with pytest.raises(ValueError, match="not valid JSON"):
        load_inputs(config=str(path), backend="ersilia")


def test_inline_annotator_errors_say_which_entry(tmp_path):
    bad = {**INLINE["annotators"][1], "cutoff": "250"}
    path = _write_config(
        tmp_path / "config.json",
        **{**INLINE, "annotators": [INLINE["annotators"][0], bad]},
    )

    with pytest.raises(ValueError, match=r"annotators\[1\] \(mw\)"):
        load_inputs(config=path, backend="ersilia")


# --- model repositories and the Hub -------------------------------------------------


def test_run_sh_needs_the_repositories_folder(inputs_dir, hub):
    with pytest.raises(ValueError, match="path-to-model-repositories"):
        _load(inputs_dir, backend="run_sh")

    assert hub == []  # a cheap local check comes before the network one


def test_the_repositories_folder_must_exist(inputs_dir, tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        _load(
            inputs_dir,
            backend="run_sh",
            path_to_model_repositories=str(tmp_path / "nowhere"),
        )


def test_every_missing_model_folder_is_listed(inputs_dir, tmp_path):
    repos = tmp_path / "repos"
    (repos / "eos9taz").mkdir(parents=True)

    with pytest.raises(ValueError) as error:
        _load(inputs_dir, backend="run_sh", path_to_model_repositories=str(repos))

    assert "['eos80sj', 'eos8nyg', 'eos4zfy']" in str(error.value)


def test_run_sh_loads_when_every_folder_is_there(inputs_dir, tmp_path):
    repos = tmp_path / "repos"
    for model in ("eos9taz", "eos80sj", "eos8nyg", "eos4zfy"):
        (repos / model).mkdir(parents=True)

    inputs = _load(inputs_dir, backend="run_sh", path_to_model_repositories=str(repos))

    assert inputs.models_dir == str(repos)
    assert inputs.generators.generators[0]._hub_model.models_dir == str(repos)


def test_the_ersilia_backend_does_not_need_the_repositories_folder(inputs_dir):
    assert _load(inputs_dir, backend="ersilia").models_dir is None


def test_a_malformed_input_never_reaches_the_hub(inputs_dir, hub):
    _write_csv(inputs_dir / "generators.csv", ["generator_id"], [["eos9taz"]])

    with pytest.raises(ValueError):
        _load(inputs_dir)

    assert hub == []


def test_a_hub_rejection_is_raised_as_is(inputs_dir, monkeypatch):
    def reject(**kwargs):
        raise ValueError("the Ersilia Model Hub does not accept these models")

    monkeypatch.setattr(config, "validate_models", reject)

    with pytest.raises(ValueError, match="does not accept"):
        _load(inputs_dir)


# --- annotators on their own --------------------------------------------------------


def test_load_annotators_reads_nothing_else(inputs_dir, hub):
    (inputs_dir / "similarity.csv").unlink()
    (inputs_dir / "generators.csv").unlink()

    specs = load_annotators(path_to_models=str(inputs_dir), backend="ersilia")

    assert [s.annotator_id for s in specs] == ["eos4zfy", "mw", "qed"]
    assert hub == [{"generators": [], "similarity": [], "annotators": ["eos4zfy"]}]


def test_in_house_annotators_need_no_repositories_even_with_run_sh(inputs_dir):
    rows = [["rdkit", "mw", "250", "higher", "primary", ""]]
    _write_csv(inputs_dir / "annotators.csv", ANNOTATOR_HEADER, rows)

    specs = load_annotators(path_to_models=str(inputs_dir), backend="run_sh")

    assert [s.annotator_id for s in specs] == ["mw"]


# --- build_annotator ----------------------------------------------------------------


def test_build_annotator_dispatches_on_source():
    assert isinstance(build_annotator("rdkit", "mw"), MolecularWeightAnnotator)
    assert isinstance(build_annotator("rdkit", "qed"), QEDAnnotator)
    hub_annotator = build_annotator("ersilia", "eos4zfy", "score", backend="ersilia")
    assert isinstance(hub_annotator, HubAnnotator) and hub_annotator.column == "score"


def test_build_annotator_rejects_unknown_names():
    with pytest.raises(ValueError, match="unknown rdkit model"):
        build_annotator("rdkit", "logp")
    with pytest.raises(ValueError, match="source must be one of"):
        build_annotator("pubchem", "mw")


# --- a workspace as the config ------------------------------------------------------


def test_inside_a_workspace_its_marker_file_is_the_config(tmp_path, monkeypatch):
    create_workspace(tmp_path)
    monkeypatch.chdir(tmp_path / "rounds")  # any folder inside the workspace will do

    inputs = load_inputs(backend="ersilia")

    root = tmp_path.resolve()
    assert inputs.similarity.model_id == "eos8nyg"
    assert inputs.models_dir == str(root / "models")
    assert inputs.envs_dir == str(root / "envs")
    assert inputs.rounds_dir == str(root / "rounds")


def test_a_flag_overrides_a_workspace_folder(tmp_path, monkeypatch):
    create_workspace(tmp_path)
    monkeypatch.chdir(tmp_path)

    inputs = load_inputs(path_to_envs=str(tmp_path / "elsewhere"), backend="ersilia")

    assert inputs.envs_dir == str(tmp_path / "elsewhere")


def test_an_explicit_config_ignores_the_workspace(tmp_path, monkeypatch):
    create_workspace(tmp_path)
    monkeypatch.chdir(tmp_path)
    path = _write_config(tmp_path / "other.json", **INLINE)

    inputs = load_inputs(config=path, backend="ersilia")

    assert inputs.envs_dir is None and inputs.rounds_dir is None


# --- environments built by fetch ----------------------------------------------------


def _repositories(tmp_path, built):
    repos, envs = tmp_path / "repos", tmp_path / "envs"
    for model in ("eos9taz", "eos80sj", "eos8nyg", "eos4zfy"):
        (repos / model).mkdir(parents=True)
    for model in built:
        env_path(envs, model).mkdir(parents=True)
        write_marker(env_path(envs, model), model=model)
    return repos, envs


def test_a_missing_or_unfinished_environment_says_to_run_fetch(inputs_dir, tmp_path):
    repos, envs = _repositories(tmp_path, built=["eos9taz", "eos80sj", "eos8nyg"])
    env_path(envs, "eos4zfy").mkdir()  # started, never marked complete

    with pytest.raises(ValueError, match=r"no complete environment for \['eos4zfy'\]"):
        _load(
            inputs_dir,
            backend="run_sh",
            path_to_model_repositories=str(repos),
            path_to_envs=str(envs),
        )


def test_complete_environments_pass_the_check(inputs_dir, tmp_path):
    repos, envs = _repositories(
        tmp_path, built=["eos9taz", "eos80sj", "eos8nyg", "eos4zfy"]
    )

    inputs = _load(
        inputs_dir,
        backend="run_sh",
        path_to_model_repositories=str(repos),
        path_to_envs=str(envs),
    )

    assert inputs.envs_dir == str(envs)


def test_a_missing_clone_with_environments_in_use_says_to_run_fetch(
    inputs_dir, tmp_path
):
    repos, envs = _repositories(tmp_path, built=[])
    (repos / "eos80sj").rmdir()

    with pytest.raises(
        ValueError, match=r"no folder for \['eos80sj'\].*chemsampler fetch"
    ):
        _load(
            inputs_dir,
            backend="run_sh",
            path_to_model_repositories=str(repos),
            path_to_envs=str(envs),
        )


# --- what fetch needs ---------------------------------------------------------------


def test_plan_fetch_lists_each_hub_model_once_with_the_workspace_folders(
    tmp_path, monkeypatch, hub
):
    create_workspace(tmp_path)
    monkeypatch.chdir(tmp_path)
    _write_csv(
        tmp_path / "input" / "generators.csv",
        ["generators"],
        [["eos9taz"], ["eos8nyg"]],
    )

    plan = plan_fetch()

    root = tmp_path.resolve()
    assert plan.model_ids == ("eos9taz", "eos8nyg", "eos4zfy")  # eos8nyg only once
    assert plan.models_dir == root / "models" and plan.envs_dir == root / "envs"
    assert plan.logs_dir == root / "logs" / "fetch"
    assert hub[0]["annotators"] == ["eos4zfy"]


def test_plan_fetch_needs_no_folder_to_exist(inputs_dir, tmp_path):
    plan = plan_fetch(
        path_to_models=str(inputs_dir),
        path_to_model_repositories=str(tmp_path / "new" / "models"),
        path_to_envs=str(tmp_path / "new" / "envs"),
    )

    assert plan.model_ids == ("eos9taz", "eos80sj", "eos8nyg", "eos4zfy")


def test_plan_fetch_outside_a_workspace_needs_the_two_folders(
    inputs_dir, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError, match="chemsampler setup"):
        plan_fetch(path_to_models=str(inputs_dir))
