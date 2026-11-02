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
    help="Also pull clones that are behind and rebuild environments that differ "
    "from their recipe.",
)
def fetch_cmd(path: str | None, update: bool) -> None:
    """Clone the project's models, download their checkpoints and build their environments."""
    from ...config import plan_fetch
    from ...hub.fetch import check_prerequisites, fetch_models, summarize
    from ...project import PROJECT_FILE, locate_project

    try:
        plan = plan_fetch(str(locate_project(path) / PROJECT_FILE))
        check_prerequisites()
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    also = ", and the rdkit environment" if plan.rdkit else ""
    click.echo(
        f"Fetching {len(plan.model_ids)} model(s): {', '.join(plan.model_ids)}{also}"
    )
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
