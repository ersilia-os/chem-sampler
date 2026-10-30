import os
import shlex
import shutil
import sys
from pathlib import Path

import click

from ..create_cli import chemsampler_cli

#: The folder made when `--path` is not given: a dummy project to try things in.
DEFAULT_NAME = "chemsampler-test"


def _folder_to_enter(root: Path) -> str:
    """`root` as one would type it after `cd`: relative when it is under here."""
    try:
        return shlex.quote(str(root.relative_to(Path.cwd())))
    except ValueError:
        return shlex.quote(str(root))


def _absolute(path) -> Path:
    return Path(os.path.abspath(os.path.expanduser(path)))


def _model_folders(models: str | None, envs: str | None) -> tuple[Path, Path]:
    """The folders to record: the ones given, and this machine's store for the rest."""
    if models is None or envs is None:
        from ... import store

        shared = store.resolve()  # a broken store config raises ValueError
        models = shared.models if models is None else models
        envs = shared.envs if envs is None else envs
    return _absolute(models), _absolute(envs)


@chemsampler_cli.command(
    "init", short_help="Create a project for one analysis (inputs and results)."
)
@click.option(
    "-p",
    "--path",
    type=click.Path(),
    help=f"Folder to create (default: ./{DEFAULT_NAME}, a dummy project).",
)
@click.option(
    "--models",
    type=click.Path(file_okay=False),
    help="Where `chemsampler fetch` clones the models. Default: this machine's "
    "store, shared by every project.",
)
@click.option(
    "--envs",
    type=click.Path(file_okay=False),
    help="Where `chemsampler fetch` builds their conda environments; prefer a "
    "local disk. Default: this machine's store.",
)
def init_cmd(path: str | None, models: str | None, envs: str | None) -> None:
    """Create a NEW project folder for one analysis's inputs and results.

    Without --path it makes ./chemsampler-test, a dummy project to try things in.
    A path that already exists is never written into: init says so and stops.
    Every location goes into the project's chemsampler.json, so the other
    commands need only the project; edit that file to change one later.
    """
    from ...project import create_project, find_project

    # abspath, not resolve: a dangling symlink must count as taken, not be followed
    root = Path(os.path.abspath(path or DEFAULT_NAME))
    if os.path.lexists(root):
        click.secho(
            f"{root} already exists; init only creates new folders. "
            "Choose another name with --path.",
            fg="red",
            err=True,
        )
        sys.exit(1)

    try:
        models_dir, envs_dir = _model_folders(models, envs)
    except ValueError as e:
        click.secho(f"could not create the project: {e}", fg="red", err=True)
        sys.exit(1)

    try:
        root.mkdir(parents=True)
        try:
            created = create_project(root, models_dir, envs_dir)
        except (OSError, ValueError):
            # the folder is ours and half-made: leave nothing that blocks a retry
            shutil.rmtree(root, ignore_errors=True)
            raise
    except (OSError, ValueError) as e:
        click.secho(f"could not create the project in {root}: {e}", fg="red", err=True)
        sys.exit(1)

    click.secho(f"Project ready in {root}:", fg="green")
    for item in created:
        click.echo(f"  {item}")
    click.echo(f"Models go to:       {models_dir}")
    click.echo(f"Environments go to: {envs_dir}")
    outer = find_project(root.parent)
    if outer is not None:
        click.secho(
            f"Note: this is inside the project {outer}; commands run inside the new "
            "folder use the new project.",
            fg="yellow",
        )
    if path is None:
        click.echo(
            f"{DEFAULT_NAME} is a dummy project; for your own, use "
            "`chemsampler init --path NAME`."
        )
    click.echo(
        f"Next: cd {_folder_to_enter(root)}, edit input/*.csv, then run "
        "`chemsampler fetch`."
    )
