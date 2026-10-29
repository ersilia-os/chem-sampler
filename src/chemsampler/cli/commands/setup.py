import sys
from pathlib import Path

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("setup")
@click.argument("path", default=".", type=click.Path(file_okay=False))
def setup_cmd(path: str) -> None:
    """Create a workspace in PATH (default: here) for inputs, models, envs and results."""
    from ...workspace import create_workspace

    root = Path(path).resolve()
    try:
        created = create_workspace(root)
    except (OSError, ValueError) as e:
        click.secho(f"could not set up {root}: {e}", fg="red")
        sys.exit(1)

    if created:
        click.secho(f"Workspace ready in {root}:", fg="green")
        for item in created:
            click.echo(f"  {item}")
    else:
        click.secho(f"{root} is already a workspace; nothing to do.", fg="green")
    click.echo(
        "Next: edit input/*.csv, then run `chemsampler fetch` from this folder. "
        "Do not move or rename it afterwards: the conda envs inside are not relocatable."
    )
