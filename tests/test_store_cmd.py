import json

import pytest
from click.testing import CliRunner

from chemsampler.cli import cli


def _store(**env):
    return CliRunner().invoke(cli, ["store"], env=env)


def test_the_store_shows_what_it_holds_and_where_it_comes_from(store_root, fake_fetch):
    fake_fetch("eos4zfy", "eos9p57")
    (store_root / "envs" / "half-built").mkdir()  # no marker: not an environment
    (store_root / "hpc.json").write_text(json.dumps({"lab": {}, "cloud": {}}))

    result = _store()

    assert result.exit_code == 0, result.output
    assert f"Store: {store_root}   (set by $CHEMSAMPLER_HOME)" in result.output
    assert "models/    2: eos4zfy, eos9p57" in result.output
    assert "envs/      2: eos4zfy, eos9p57" in result.output
    assert "hpc.json   profiles: cloud, lab" in result.output


def test_the_store_says_how_to_move_it(store_root):
    result = _store()

    assert "export CHEMSAMPLER_HOME=/path/to/.chemsampler" in result.output
    assert ">> ~/.bashrc" in result.output
    assert "environments cannot be" in result.output


def test_without_the_variable_the_default_is_shown_and_not_made(tmp_path):
    result = _store(CHEMSAMPLER_HOME="", HOME=str(tmp_path))

    assert result.exit_code == 0, result.output
    expected = f"Store: {tmp_path / '.chemsampler'}   (the default: $CHEMSAMPLER_HOME"
    assert expected in result.output
    assert "not created yet" in result.output
    assert not (tmp_path / ".chemsampler").exists()


def test_a_malformed_profiles_file_is_reported(store_root):
    store_root.mkdir()
    (store_root / "hpc.json").write_text("[1]")

    result = _store()

    assert result.exit_code == 0, result.output
    assert "profiles: malformed:" in result.output


def test_an_old_config_file_is_refused(store_root):
    store_root.mkdir()
    (store_root / "config.json").write_text("{}")

    result = _store()

    assert result.exit_code == 1
    assert "no longer read" in result.output


@pytest.mark.parametrize("name", sorted(cli.commands))
def test_each_command_is_described_in_one_paragraph(name):
    # details belong in the README and docs/, not in --help
    assert "\n\n" not in (cli.commands[name].help or "").strip()
