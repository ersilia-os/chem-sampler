import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("fetch")
@click.option(
    "--config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="JSON file with the paths and the model lists. Inside a project "
    "(see `chemsampler init`) its chemsampler.json is used by default.",
)
@click.option(
    "--path-to-models",
    default=None,
    type=click.Path(file_okay=False),
    help="Folder with similarity.csv, generators.csv and annotators.csv.",
)
@click.option(
    "--path-to-model-repositories",
    default=None,
    type=click.Path(file_okay=False),
    help="Where the models are cloned, one folder per model. Default: this "
    "machine's store (~/.chemsampler/models, see `chemsampler catalog`).",
)
@click.option(
    "--path-to-envs",
    default=None,
    type=click.Path(file_okay=False),
    help="Where the conda environments are built, one per model. Default: this "
    "machine's store (~/.chemsampler/envs). Keep it on a local disk.",
)
@click.option(
    "--similarity",
    "similarity_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with a 'similarity' column. Overrides the config file and --path-to-models.",
)
@click.option(
    "--generators",
    "generators_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with a 'generators' column. Overrides the config file and --path-to-models.",
)
@click.option(
    "--annotators",
    "annotators_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with the annotators. Overrides the config file and --path-to-models.",
)
@click.option(
    "--update",
    is_flag=True,
    default=False,
    help="Pull clones that are behind (fast-forward only). A clone with local "
    "changes is never touched. Without it, existing clones are only inspected.",
)
@click.option(
    "--rebuild",
    is_flag=True,
    default=False,
    help="Rebuild the environments that differ from their model's recipe. "
    "Without it, differences are only reported.",
)
def fetch_cmd(
    config: str | None,
    path_to_models: str | None,
    path_to_model_repositories: str | None,
    path_to_envs: str | None,
    similarity_file: str | None,
    generators_file: str | None,
    annotators_file: str | None,
    update: bool,
    rebuild: bool,
) -> None:
    """Clone the models, download their checkpoints and build their conda environments."""
    from ...config import plan_fetch
    from ...hub.fetch import check_prerequisites, fetch_models, summarize

    try:
        plan = plan_fetch(
            config,
            path_to_models,
            path_to_model_repositories,
            path_to_envs,
            similarity=similarity_file,
            generators=generators_file,
            annotators=annotators_file,
        )
        check_prerequisites()
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"Fetching {len(plan.model_ids)} model(s): {', '.join(plan.model_ids)}")
    click.echo(f"  clones:       {plan.models_dir}")
    click.echo(f"  environments: {plan.envs_dir}")
    click.echo(f"  logs:         {plan.logs_dir}")
    results = fetch_models(plan, update=update, rebuild=rebuild, say=click.echo)

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
