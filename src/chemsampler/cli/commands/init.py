import sys
from pathlib import Path

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command(
    "init", short_help="Create a project for one analysis (inputs and results)."
)
@click.argument("path", default=".", type=click.Path(file_okay=False))
def init_cmd(path: str) -> None:
    """Create a project in PATH (default: here) for one analysis's inputs and results."""
    from ...project import create_project

    root = Path(path).resolve()
    try:
        created = create_project(root)
    except (OSError, ValueError) as e:
        click.secho(f"could not create the project in {root}: {e}", fg="red")
        sys.exit(1)

    if created:
        click.secho(f"Project ready in {root}:", fg="green")
        for item in created:
            click.echo(f"  {item}")
    else:
        click.secho(f"{root} is already a project; nothing to do.", fg="green")
    click.echo(
        "Next: edit input/*.csv, then run `chemsampler fetch` from this folder. "
        "The models themselves go to this machine's store, shared by every project "
        "(`chemsampler catalog` shows it)."
    )
