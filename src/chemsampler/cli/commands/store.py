import os
import sys

import click

from ..create_cli import chemsampler_cli

_PLACEHOLDER = "/path/to/.chemsampler"


@chemsampler_cli.command(
    "store", short_help="Show where the models and environments are kept."
)
def store_cmd() -> None:
    """Show the store (models, environments, cluster profiles) and how to move it."""
    from ... import store
    from ...hpc.profiles import HPC_FILE
    from ...hub.envs import is_complete
    from ...utils.jsonfile import read_object

    try:
        shared = store.resolve()
    except ValueError as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    variable = store.HOME_VARIABLE
    source = (
        f"set by ${variable}"
        if os.environ.get(variable)
        else f"the default: ${variable} is not set"
    )
    click.echo(f"Store: {shared.root}   ({source})")
    if not shared.root.is_dir():
        click.echo("  not created yet; `chemsampler fetch` makes it")
    else:
        models = _folders(shared.models)
        envs = [p.name for p in _paths(shared.envs) if is_complete(p)]
        click.echo(f"  models/    {_listing(models)}")
        click.echo(f"  envs/      {_listing(envs)}")
        hpc = shared.root / HPC_FILE
        if not hpc.is_file():
            profiles = "none"
        else:
            try:
                profiles = ", ".join(sorted(read_object(hpc))) or "none"
            except ValueError as e:
                profiles = f"malformed: {e}"
        click.echo(f"  {HPC_FILE}   profiles: {profiles}")

    click.echo()
    export = f"export {variable}={_PLACEHOLDER}"
    keep = f"echo '{export}' >> ~/.bashrc"
    click.echo("To use another folder:")
    click.echo(f"  {export:<{len(keep)}}   # this terminal")
    click.echo(f"  {keep}   # every new terminal")
    click.echo("Then run `chemsampler fetch` in your projects: environments cannot be")
    click.echo("moved, they record their own path. Projects made before keep the")
    click.echo("folders their chemsampler.json names.")


def _paths(folder):
    if not folder.is_dir():
        return []
    return sorted(p for p in folder.iterdir() if p.is_dir())


def _folders(folder) -> list[str]:
    return [p.name for p in _paths(folder)]


def _listing(names: list[str]) -> str:
    return f"{len(names)}: {', '.join(names)}" if names else "none"
