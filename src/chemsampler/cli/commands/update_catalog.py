import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("update-catalog")
def update_catalog_cmd() -> None:
    """Rebuild the model catalog from the Ersilia Model Hub (about 20 s)."""
    from ...hub import catalog

    click.echo("Querying the Ersilia Model Hub (about 20 s)...")
    try:
        result = catalog.update_catalog()
    except RuntimeError as e:
        click.secho(f"{e}\nThe existing catalog was left untouched.", fg="red")
        sys.exit(1)
    except OSError as e:
        click.secho(
            f"could not write {catalog.CATALOG_PATH}: {e}\n"
            "Run update-catalog from a development checkout (pip install -e .).",
            fg="red",
        )
        sys.exit(1)

    counts = ", ".join(f"{result.catalog.count(r)} {r}" for r in catalog.ROLES)
    if not result.written:
        click.secho(
            f"Catalog is up to date: {counts} (built {result.catalog.built}).",
            fg="green",
        )
        return

    if result.previous is None:
        click.secho(f"Catalog created: {counts}.", fg="green")
    else:
        click.secho(
            f"Catalog updated (was built {result.previous.built}): {counts}.",
            fg="green",
        )
        for label, entries in (
            ("new", result.diff.added),
            ("removed", result.diff.removed),
        ):
            for role in catalog.ROLES:
                models = [e for e in entries if e.role == role]
                if models:
                    listed = ", ".join(f"{e.identifier} ({e.title})" for e in models)
                    click.echo(f"  {label} {role}: {listed}")
    click.echo(f"Commit {catalog.CATALOG_PATH} to share it.")
