import sys

import click

from ..create_cli import chemsampler_cli
from ..options import confirm, project_option, split_ids

#: The id of the rdkit environment, as `fetch` names it.
_RDKIT = "rdkit"


@chemsampler_cli.command(
    "remove", short_help="Delete models (clone and environment) from the store."
)
@project_option
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Do not ask, and remove clones that hold local changes or commits.",
)
@click.option(
    "--unused",
    is_flag=True,
    default=False,
    help="Remove every model the project's input lists do not name.",
)
@click.argument("ids", nargs=-1)
def remove_cmd(
    path: str | None, force: bool, unused: bool, ids: tuple[str, ...]
) -> None:
    """Remove these models (or `rdkit`) from the store; the input lists stay as they are."""
    from ...config import plan_fetch, store_folders
    from ...hub.remove import plan_removal, plan_unused, remove
    from ...models.annotator import RDKIT_ANNOTATORS
    from ...project import PROJECT_FILE, locate_project

    names = split_ids(ids)
    if bool(names) == unused:
        raise click.UsageError("give model ids, or --unused: one of the two")
    in_house = [n for n in names if n in RDKIT_ANNOTATORS]
    if in_house:
        click.secho(
            f"{', '.join(in_house)}: not models but RDKit annotators, which run in "
            f"the rdkit environment: use `chemsampler remove {_RDKIT}`",
            fg="red",
        )
        sys.exit(1)

    try:
        config = str(locate_project(path) / PROJECT_FILE)
        models_dir, envs_dir, _ = store_folders(config)
    except (ValueError, RuntimeError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    used, rdkit_used = (), False
    try:
        inputs = plan_fetch(config)
        used, rdkit_used = inputs.model_ids, bool(inputs.rdkit_annotators)
    except ValueError as e:
        if unused:  # what is unused cannot be judged without the lists
            click.secho(f"input/*.csv could not be read: {e}", fg="red")
            sys.exit(1)
        click.secho(f"(this project's inputs could not be read: {e})", fg="yellow")

    if unused:
        plan = plan_unused(models_dir, envs_dir, used, rdkit_used, force)
    else:
        plan = plan_removal(models_dir, envs_dir, names, force)

    for model_id in plan.absent:
        click.echo(f"{model_id}: not in the store")
    for item in plan.blocked:
        click.secho(
            f"{item.model_id}: {item.blocked}, not removed (use -f to remove anyway)",
            fg="yellow",
        )
    if not plan.to_remove:
        if plan.removals or plan.absent:
            click.echo("Nothing removed.")
            sys.exit(1)
        click.echo("Nothing unused.")
        return

    click.echo(f"Removing from {models_dir} and {envs_dir}:")
    width = max(len(r.model_id) for r in plan.to_remove)
    for item in plan.to_remove:
        click.echo(f"  {item.model_id:<{width}}  {item.parts}")
    in_use = [
        r.model_id for r in plan.to_remove if _in_use(r.model_id, used, rdkit_used)
    ]
    for model_id in in_use:
        click.echo(
            f"{model_id} is used by this project: `chemsampler run` will need "
            "`chemsampler fetch` again."
        )
    click.echo(
        "Other projects sharing this store may use these models; chemsampler "
        "cannot tell."
    )
    if not force:
        try:
            agreed = confirm(
                "Remove these?",
                "remove would delete everything listed above; run it with -f to "
                "remove without asking",
            )
        except ValueError as e:
            click.secho(str(e), fg="red")
            sys.exit(1)
        if not agreed:
            click.echo("Nothing removed.")
            return

    failures = remove(plan, say=click.echo)
    done = [r.model_id for r in plan.to_remove if r.model_id not in dict(failures)]
    if done:
        click.echo(f"Removed {len(done)} model(s): {', '.join(done)}")
    for model_id, reason in failures:
        click.secho(f"{model_id}: could not be removed: {reason}", fg="red")
    if failures or plan.blocked:
        sys.exit(1)


def _in_use(model_id: str, used, rdkit_used: bool) -> bool:
    return rdkit_used if model_id == _RDKIT else model_id in used
