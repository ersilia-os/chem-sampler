import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option


@chemsampler_cli.command(
    "resume", short_help="Continue a stopped run from the rounds it finished."
)
@project_option
@click.argument("run", required=False)
def resume_cmd(path: str | None, run: str | None) -> None:
    """Continue a run that stopped: Ctrl-C, a failed model, a closed terminal, a
    preempted or cancelled cluster job.

    RUN is one of the project's runs (a folder of rounds/, or its path); without
    it, the latest. The rounds it finished are read back from its folder, not run
    again, and the one it was in is run again whole, with the inputs and options
    the run started with (including --hpc).
    """
    from ... import runs
    from ...config import load_inputs
    from ...project import PROJECT_FILE, Project, locate_project
    from ...utils.logging import logger
    from .run import execute

    try:
        root = locate_project(path)
        run_dir = runs.find_run(Project.load(root).rounds, run)
        record = runs.read_record(run_dir)
        if not (run_dir / runs.INPUT_FOLDER).is_dir() or "status" not in record.get(
            "run", {}
        ):
            raise ValueError(
                f"{run_dir} was made before chemsampler could resume runs: it has "
                "no copy of its inputs"
            )
        if record["run"]["status"] == "finished":
            click.secho(f"{run_dir.name} already finished: see {run_dir}", fg="green")
            return
        if not record["run"]["verbose"]:
            logger.set_quiet_mode(True)
        # the inputs as the run found them, not as the project has them now
        inputs = load_inputs(
            str(root / PROJECT_FILE),
            path_to_models=str(run_dir / runs.INPUT_FOLDER),
            path_to_model_repositories=record["path_to_model_repositories"],
            path_to_envs=record["path_to_envs"],
        )
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"Resuming {run_dir.name} ({record['run']['status']}).")
    execute(inputs, record, run_dir)
