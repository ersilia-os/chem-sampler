import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option


@chemsampler_cli.command(
    "fetch", short_help="Clone models, download checkpoints, build their conda envs."
)
@project_option
@click.option(
    "--update",
    is_flag=True,
    default=False,
    help="Bring existing models up to date: pull clones that are behind "
    "(fast-forward only; a clone with local changes is never touched) and rebuild "
    "environments that differ from their recipe. Without it, they are only checked.",
)
def fetch_cmd(path: str | None, update: bool) -> None:
    """Get the models a project needs: clone each one, download its checkpoints
    and build its conda environment, in the folders the project names.
    """
    from ...config import plan_fetch
    from ...hub.fetch import check_prerequisites, fetch_models, summarize
    from ...project import PROJECT_FILE, locate_project

    try:
        plan = plan_fetch(str(locate_project(path) / PROJECT_FILE))
        check_prerequisites()
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"Fetching {len(plan.model_ids)} model(s): {', '.join(plan.model_ids)}")
    click.echo(f"  clones:       {plan.models_dir}")
    click.echo(f"  environments: {plan.envs_dir}")
    click.echo(f"  logs:         {plan.logs_dir}")
    results = fetch_models(plan, update=update, rebuild=update, say=click.echo)

    click.echo()
    for line in summarize(results):
        click.echo(line)
    failed = [r.model_id for r in results if not r.ok]
    if failed:
        click.secho(f"{len(failed)} model(s) failed: {', '.join(failed)}", fg="red")
        sys.exit(1)
    differing = [r.model_id for r in results if r.env == "differs"]
    if differing:
        click.secho(
            "Done, with environments that differ from their recipe.", fg="yellow"
        )
    else:
        click.secho("All models are ready.", fg="green")
