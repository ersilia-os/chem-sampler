import contextlib
import logging
import sys

import click

from ..create_cli import chemsampler_cli

_ROLES = ["generators", "similarity", "annotators"]


@chemsampler_cli.command("catalog")
@click.option(
    "--hub",
    is_flag=True,
    default=False,
    help="List what the Ersilia Model Hub offers (from the bundled catalog) "
    "instead of what this machine has.",
)
@click.option(
    "--role",
    type=click.Choice(_ROLES),
    default=None,
    help="Only this kind of model.",
)
@click.option(
    "--update",
    is_flag=True,
    default=False,
    help="Rebuild the bundled catalog from the Hub (about 20 s, needs a network) "
    "and report which models are new or gone.",
)
@click.option(
    "--config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="JSON file with the paths. Inside a project its chemsampler.json is "
    "used by default.",
)
@click.option(
    "--path-to-model-repositories",
    default=None,
    type=click.Path(file_okay=False),
    help="The folder with the clones. Default: this machine's store.",
)
@click.option(
    "--path-to-envs",
    default=None,
    type=click.Path(file_okay=False),
    help="The folder with the environments. Default: this machine's store.",
)
def catalog_cmd(
    hub: bool,
    role: str | None,
    update: bool,
    config: str | None,
    path_to_model_repositories: str | None,
    path_to_envs: str | None,
) -> None:
    """List the models on this machine, or those the Hub offers, or refresh the catalog."""
    if update:
        if hub or role:
            click.secho("--update cannot be combined with --hub or --role.", fg="red")
            sys.exit(2)
        _update()
        return

    from ...config import locate_models

    try:
        models_dir, envs_dir = locate_models(
            config, path_to_model_repositories, path_to_envs
        )
    except ValueError as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    if hub:
        _show_hub(models_dir, role)
    else:
        _show_local(
            models_dir, envs_dir, role, config, path_to_model_repositories, path_to_envs
        )


@contextlib.contextmanager
def _quiet():
    """Hide the progress lines of the loaders: a listing should print only the listing."""
    from ...utils.logging import logger

    previous = logger.logger.level
    logger.logger.setLevel(logging.WARNING)
    try:
        yield
    finally:
        logger.logger.setLevel(previous)


def _read_catalog():
    """The bundled catalog, or None (with a hint) if it cannot be read."""
    from ...hub import catalog as hub_catalog

    try:
        return hub_catalog.read_catalog()
    except ValueError as e:
        click.secho(str(e), fg="yellow")
        return None


def _show_local(models_dir, envs_dir, role, config, repos_flag, envs_flag) -> None:
    from ...config import plan_fetch
    from ...hub.listing import list_local
    from ...project import find_project

    needed = ()
    try:
        with _quiet():
            needed = plan_fetch(config, None, repos_flag, envs_flag).model_ids
    except ValueError as e:
        # outside a project there are no inputs to read, which is fine; inside
        # one, the user should know why nothing is marked
        if config is not None or find_project() is not None:
            click.secho(f"(this project's inputs could not be read: {e})", fg="yellow")

    models = list_local(
        models_dir,
        envs_dir,
        needed=needed,
        catalog=_read_catalog(),
        role=role,
    )
    where = f"Models in {models_dir}"
    if envs_dir is not None:
        where += f" (environments in {envs_dir})"
    click.echo(where)

    if not models:
        click.echo(
            "No models here yet. Run `chemsampler fetch` in a project to get some, "
            "or `chemsampler catalog --hub` to see what exists."
        )
        return

    width = max(len(m.model_id) for m in models)
    click.echo(f"  {'role':<11} {'model':<{width}}  {'env':<10} title")
    for m in models:
        mark = "*" if m.needed else " "
        click.echo(
            f"{mark} {m.role:<11} {m.model_id:<{width}}  {m.env:<10} {m.title}".rstrip()
        )
    if any(m.needed for m in models):
        click.echo("* needed by this project")
    lacking = [m.model_id for m in models if m.needed and not m.cloned]
    unbuilt = [
        m.model_id
        for m in models
        if m.needed and m.cloned and m.env in ("missing", "unfinished")
    ]
    if lacking or unbuilt:
        click.secho(
            f"{len(lacking) + len(unbuilt)} model(s) this project needs are not ready "
            f"({', '.join(lacking + unbuilt)}); run `chemsampler fetch`.",
            fg="yellow",
        )


def _show_hub(models_dir, role) -> None:
    from ...hub.catalog import ROLES
    from ...hub.listing import list_hub

    catalog = _read_catalog()
    if catalog is None:
        sys.exit(1)
    entries = list_hub(catalog, models_dir, role)
    click.echo(
        f"The Ersilia Model Hub, as of the catalog built {catalog.built} "
        f"(ersilia {catalog.ersilia_version})"
    )
    if role is None:
        for each in ROLES:
            present = sum(1 for e in entries if e.role == each and e.in_store)
            click.echo(
                f"  {each:<11} {catalog.count(each):>4} models, {present} on this machine"
            )
        click.echo("List one kind with --role generators|similarity|annotators.")
        return

    width = max((len(e.model_id) for e in entries), default=5)
    for e in entries:
        notes = []
        if e.in_store:
            notes.append("on this machine")
        if not e.supported:
            notes.append("not supported by chemsampler")
        suffix = f"  [{', '.join(notes)}]" if notes else ""
        click.echo(f"  {e.model_id:<{width}}  {e.title}{suffix}")
    click.echo(f"{len(entries)} {role}.")


def _update() -> None:
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
            "Run `catalog --update` from a development checkout (pip install -e .).",
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
            for each in catalog.ROLES:
                models = [e for e in entries if e.role == each]
                if models:
                    listed = ", ".join(f"{e.identifier} ({e.title})" for e in models)
                    click.echo(f"  {label} {each}: {listed}")
    click.echo(f"Commit {catalog.CATALOG_PATH} to share it.")
