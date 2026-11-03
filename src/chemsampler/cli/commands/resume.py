import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option

#: Options of older runs that a run can no longer take, so cannot resume with.
_REMOVED_OPTIONS = ("tanimoto_cutoff", "original_seed_smiles")


@chemsampler_cli.command(
    "resume", short_help="Continue a stopped run from the rounds it finished."
)
@project_option
@click.argument("run", required=False)
@click.option(
    "--detach",
    is_flag=True,
    default=False,
    help="Continue it as a cluster job (runs that used --hpc).",
)
def resume_cmd(path: str | None, run: str | None, detach: bool) -> None:
    """Continue a stopped run from the rounds it finished (default: the latest run)."""
    from ... import runs
    from ...config import load_inputs
    from ...project import PROJECT_FILE, Project, locate_project
    from ...utils.logging import logger
    from .run import detach_run, execute

    try:
        root = locate_project(path)
        run_dir = runs.find_run(Project.load(root).runs, run)
        record = runs.read_record(run_dir)
        if not (run_dir / runs.INPUT_FOLDER).is_dir() or "status" not in record.get(
            "run", {}
        ):
            raise ValueError(
                f"{run_dir} was made before chemsampler could resume runs: it has "
                "no copy of its inputs"
            )
        if any(record["run"].get(k) is not None for k in _REMOVED_OPTIONS):
            raise ValueError(
                f"{run_dir.name} was made by an older chemsampler with "
                "--tanimoto-cutoff or --original-seed-smiles, which runs no longer "
                "take (a Tanimoto limit is now an annotator, rdkit,tanimoto): start "
                "it again"
            )
        if record["run"]["status"] == "finished":
            click.secho(f"{run_dir.name} already finished: see {run_dir}", fg="green")
            return
        if "seeds" not in record["run"]:
            raise ValueError(
                f"{run_dir.name} was made by an older chemsampler, before each "
                "seed had its own folder: start it again"
            )
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

    if detach:
        detach_run(inputs, record, run_dir, root)
        return
    click.echo(f"Resuming {run_dir.name} ({record['run']['status']}).")
    execute(inputs, record, run_dir)
