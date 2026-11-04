import os
import signal
import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option


def _chemsampler_version() -> str:
    import importlib.metadata

    try:
        return importlib.metadata.version("chemsampler")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _stop_on_sigterm(signum, frame):
    # Slurm's preemption and scancel: stop as Ctrl-C does, cancelling the jobs
    raise KeyboardInterrupt


def execute(inputs, record: dict, run_dir, input_dir=None) -> None:
    """
    Run, or continue, the run in `run_dir`, then report and exit as `run` does.

    Each seed is climbed in turn, in its own folder (`runs.seed_folders`). Every
    round's table is saved as soon as it is made, and rounds already saved are
    read back instead of run again (see `hill_climb`'s `checkpoint`), so a seed
    that finished is only read. A stop or a failure saves the current seed's
    finished rounds and leaves the later seeds for `resume`. `input_dir` is given
    for a new run: the folder is made, with a copy of the inputs, once the
    cluster (if any) has been checked.
    """
    import contextlib
    import time

    import pandas as pd

    from ... import execution, runs
    from ...optimize import hill_climb, write_results

    options = record["run"]
    hpc = options["hpc"]
    name = os.path.basename(run_dir)
    start_time = time.time()
    try:
        cluster, job_alive = contextlib.nullcontext(), None
        if hpc is not None:
            from ...hpc import prepare

            cluster = prepare(hpc, inputs, run_dir)  # checks first, submits nothing
            job_alive = cluster.job_alive
            click.echo(
                f"Running every model on the cluster of profile {hpc!r}; jobs are "
                "checked every 10 s (every 60 s once a step passes 5 min)."
            )
        if input_dir is not None:
            runs.start(run_dir, input_dir, record)
        if job_alive is not None:
            runs.check_driver(record, job_alive)
        runs.take_lock(run_dir, job_alive)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)
    except KeyboardInterrupt:
        click.secho("\nInterrupted before the first round.", fg="yellow")
        sys.exit(0)

    runs.set_status(run_dir, record, "running")
    seeds = options["seeds"]
    folders = runs.seed_folders(run_dir, len(seeds))
    summaries = []
    folder = folders[0]
    previous = signal.signal(signal.SIGTERM, _stop_on_sigterm)
    try:
        for k, (seed, folder) in enumerate(zip(seeds, folders), start=1):
            if len(seeds) > 1:
                click.echo(f"Seed {k} of {len(seeds)}: {seed} ({folder.name}/)")
            # each seed's jobs and jobs.json go in its own folder
            here = cluster if hpc is None else cluster.for_folder(folder)
            with here as executor, execution.use(executor):
                summary, candidates_by_round = hill_climb(
                    inputs.generators,
                    inputs.annotators,
                    similarity=inputs.similarity,
                    seed_smiles=seed,
                    n_rounds=options["n_rounds"],
                    tolerance=options["tolerance"],
                    checkpoint=folder,
                )
            write_results(summary, candidates_by_round, folder)
            summaries.append(summary)
        elapsed = time.time() - start_time
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        # a model that failed midway: the rounds before it are still results
        if hasattr(e, "partial"):
            write_results(*e.partial, folder)
            click.secho(f"Finished rounds written to {folder}", fg="yellow")
        runs.set_status(run_dir, record, "failed")
        runs.release_lock(run_dir)
        click.echo(f"Once fixed, continue it with `chemsampler resume {name}`.")
        sys.exit(1)
    except KeyboardInterrupt as e:
        click.secho("\nInterrupted! Saving partial results...", fg="yellow")
        write_results(*getattr(e, "partial", (pd.DataFrame(), {})), folder)
        runs.set_status(run_dir, record, "stopped")
        runs.release_lock(run_dir)
        click.secho(f"Partial results written to {folder}", fg="yellow")
        click.echo(f"Continue it with `chemsampler resume {name}`.")
        sys.exit(0)
    finally:
        signal.signal(signal.SIGTERM, previous)

    runs.set_status(run_dir, record, "finished")
    runs.release_lock(run_dir)
    _report(seeds, summaries, folders, elapsed, options["verbose"], run_dir)


def detach_run(inputs, record: dict, run_dir, project_root, input_dir=None) -> None:
    """
    Send the run in `run_dir` to the cluster as a job of its own, and say so.

    `input_dir` is given for a new run, whose folder is made once the cluster has
    been checked. A run that something still runs, or that a job of its own
    waits to run, is refused.
    """
    from ... import runs
    from ...hpc import say_detached, submit_driver

    def ready(executor):
        if input_dir is not None:
            record["run"]["status"] = "queued"
            runs.start(run_dir, input_dir, record)
            return
        runs.check_driver(record, executor.job_alive)
        runs.check_lock(run_dir, executor.job_alive)
        runs.release_lock(run_dir)  # left by a run that is gone

    hpc = record["run"]["hpc"]
    if hpc is None:
        click.secho(
            "--detach sends the run to a cluster as a job: it needs --hpc NAME",
            fg="red",
        )
        sys.exit(1)
    try:
        job_id = submit_driver(hpc, inputs, run_dir, project_root, before_submit=ready)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        if input_dir is not None and os.path.isdir(run_dir):
            import shutil

            shutil.rmtree(run_dir)  # a new run that never left: nothing to resume
        sys.exit(1)
    record["run"]["driver_job"] = job_id
    runs.set_status(run_dir, record, "queued")
    say_detached(run_dir, job_id)


def _report(seeds, summaries, folders, elapsed: float, verbose: bool, run_dir) -> None:
    for seed, summary, folder in zip(seeds, summaries, folders):
        where = f" from {seed} ({folder.name}/)" if len(seeds) > 1 else ""
        if summary.empty:
            click.secho(f"No candidates were generated{where}.", fg="yellow")
            continue
        best = summary.iloc[-1]
        click.secho(
            f"Best candidate{where}: {best['smiles']} (round {best['round']}, "
            f"score={best['score']})",
            fg="green",
        )
        if not verbose:
            top = summary.nlargest(min(5, len(summary)), "score")
            click.secho("Top scores:", fg="cyan")
            for _, row in top.iterrows():
                click.secho(
                    f"  Round {int(row['round'])}: {row['score']:.2f}", fg="cyan"
                )
    if not verbose:
        click.secho(f"Total runtime: {elapsed / 60:.1f} min", fg="cyan")
    click.secho(f"Results written to {run_dir}", fg="green")


@chemsampler_cli.command(
    "run", short_help="Generate and score candidates around a seed, over rounds."
)
@project_option
@click.option(
    "--seed-smiles",
    default=None,
    help="One SMILES to start from (default: the seeds in input/smiles.csv).",
)
@click.option(
    "--n-rounds",
    default=5,
    show_default=True,
    type=int,
    help="Most rounds.",
)
@click.option(
    "--tolerance",
    default=0.0,
    show_default=True,
    type=float,
    help="How much a round must improve on the best to count.",
)
@click.option(
    "--hpc",
    metavar="NAME",
    default=None,
    help="Run the models as Slurm jobs with this cluster profile (see docs/hpc.md).",
)
@click.option(
    "--detach",
    is_flag=True,
    default=False,
    help="With --hpc: make the run itself a cluster job and return.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Show every model call, with timestamps.",
)
def run_cmd(
    path: str | None,
    seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    hpc: str | None,
    detach: bool,
    verbose: bool,
) -> None:
    """Generate and score candidates around each seed, round by round."""
    from ...config import load_inputs
    from ...project import PROJECT_FILE, Project, locate_project
    from ...utils.logging import logger

    if not verbose:
        logger.set_quiet_mode(True)

    try:
        root = locate_project(path)
        inputs = load_inputs(str(root / PROJECT_FILE))
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)
    seeds = inputs.seeds
    if seed_smiles is not None:
        from rdkit import Chem

        # RDKit would read "CCO CCN" as CCO, the rest being a name
        spaced = any(c.isspace() for c in seed_smiles.strip())
        if spaced or Chem.MolFromSmiles(seed_smiles) is None:
            click.secho(
                f"--seed-smiles takes one SMILES; RDKit cannot read {seed_smiles!r} "
                "as one (put several seeds in input/smiles.csv)",
                fg="red",
            )
            sys.exit(1)
        if seeds:
            click.secho("--seed-smiles given: ignoring input/smiles.csv", fg="yellow")
        seeds = (seed_smiles,)
    if not seeds:
        click.secho(
            "give a seed: input/smiles.csv (one SMILES per row) or --seed-smiles",
            fg="red",
        )
        sys.exit(1)
    # the models it used and how it was launched: enough to run it again
    record = {
        **inputs.describe(),
        "run": {
            "chemsampler_version": _chemsampler_version(),
            "seeds": list(seeds),
            "n_rounds": n_rounds,
            "tolerance": tolerance,
            "hpc": hpc,
            "verbose": verbose,
            "status": "running",
        },
    }
    run_dir, input_dir = inputs.new_run_dir(), Project.load(root).input
    if detach:
        detach_run(inputs, record, run_dir, root, input_dir)
    else:
        execute(inputs, record, run_dir, input_dir)
