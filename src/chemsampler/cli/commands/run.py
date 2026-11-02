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

    Every round's table is saved as soon as it is made, and rounds already in
    `run_dir` are read back instead of run again (see `hill_climb`'s
    `checkpoint`). `input_dir` is given for a new run: the folder is made, with a
    copy of the inputs, once the cluster (if any) has been checked.
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
                f"checked every {cluster.profile.poll_seconds:g} s."
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
    summary, candidates_by_round = pd.DataFrame(), {}
    previous = signal.signal(signal.SIGTERM, _stop_on_sigterm)
    try:
        with cluster as executor, execution.use(executor):
            summary, candidates_by_round = hill_climb(
                inputs.generators,
                inputs.annotators,
                mode=options["mode"],
                similarity=inputs.similarity,
                seed_smiles=options["seed_smiles"],
                original_seed_smiles=options["original_seed_smiles"],
                n_rounds=options["n_rounds"],
                n_steps=options["n_steps"],
                tolerance=options["tolerance"],
                tanimoto_cutoff=options["tanimoto_cutoff"],
                tanimoto_direction=options["tanimoto_direction"],
                checkpoint=run_dir,
            )
        elapsed = time.time() - start_time
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        # a model that failed midway: the rounds before it are still results
        if hasattr(e, "partial"):
            write_results(*e.partial, run_dir)
            click.secho(f"Finished rounds written to {run_dir}", fg="yellow")
        runs.set_status(run_dir, record, "failed")
        runs.release_lock(run_dir)
        click.echo(f"Once fixed, continue it with `chemsampler resume {name}`.")
        sys.exit(1)
    except KeyboardInterrupt as e:
        click.secho("\nInterrupted! Saving partial results...", fg="yellow")
        summary, candidates_by_round = getattr(
            e, "partial", (summary, candidates_by_round)
        )
        write_results(summary, candidates_by_round, run_dir)
        runs.set_status(run_dir, record, "stopped")
        runs.release_lock(run_dir)
        click.secho(f"Partial results written to {run_dir}", fg="yellow")
        click.echo(f"Continue it with `chemsampler resume {name}`.")
        sys.exit(0)
    finally:
        signal.signal(signal.SIGTERM, previous)

    write_results(summary, candidates_by_round, run_dir)
    runs.set_status(run_dir, record, "finished")
    runs.release_lock(run_dir)
    _report(summary, elapsed, options["verbose"], run_dir)


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


def _report(summary, elapsed: float, verbose: bool, run_dir) -> None:
    if summary.empty:
        click.secho("No candidates were generated; nothing to report.", fg="yellow")
    else:
        best = summary.iloc[-1]
        click.secho(
            f"Best candidate: {best['smiles']} (round {best['round']}, "
            f"score={best['score']})",
            fg="green",
        )

        if not verbose:
            click.echo()
            runtime_mins = elapsed / 60
            click.secho(f"Total runtime: {runtime_mins:.1f} min", fg="cyan")

            top_n = min(5, len(summary))
            top_scores = (
                summary.nlargest(top_n, "score")
                if "score" in summary.columns
                else summary.head(top_n)
            )
            if len(top_scores) > 0:
                click.secho("Top scores:", fg="cyan")
                for idx, row in top_scores.iterrows():
                    click.secho(
                        f"  Round {int(row['round'])}: {row['score']:.2f}",
                        fg="cyan",
                    )

    click.secho(f"Results written to {run_dir}", fg="green")


@chemsampler_cli.command(
    "run", short_help="Generate and score candidates around a seed, over rounds."
)
@project_option
@click.option(
    "--mode",
    type=click.Choice(["joint", "sequential", "weighted", "incremental"]),
    required=True,
    help="How the primary annotators combine into a round's winner.",
)
@click.option("--seed-smiles", default=None, help="SMILES of the starting molecule.")
@click.option(
    "--original-seed-smiles",
    default=None,
    help="SMILES of the true original molecule in a manually re-seeded chain "
    "of runs. Adds a tanimoto_to_original_seed column, independent of "
    "--seed-smiles.",
)
@click.option("--n-rounds", default=5, show_default=True, type=int)
@click.option(
    "--n-steps",
    default=None,
    type=int,
    help="Number of interpolation levels per annotator. Required for "
    "--mode incremental; unused (and rejected) otherwise.",
)
@click.option("--tolerance", default=0.0, show_default=True, type=float)
@click.option(
    "--tanimoto-cutoff",
    default=None,
    type=float,
    help="Similarity cutoff vs --seed-smiles, gated by --tanimoto-direction. "
    "Requires --seed-smiles.",
)
@click.option(
    "--tanimoto-direction",
    type=click.Choice(["higher", "lower"]),
    default="higher",
    show_default=True,
    help="'higher' keeps candidates at least this similar to the seed; "
    "'lower' pushes toward novelty (at most this similar).",
)
@click.option(
    "--hpc",
    metavar="NAME",
    default=None,
    help="Run every model as a Slurm job, with the cluster profile NAME in "
    "hpc.json in chemsampler's folder ($CHEMSAMPLER_HOME, default ~/.chemsampler). "
    "Default: here, one model at a time.",
)
@click.option(
    "--detach",
    is_flag=True,
    default=False,
    help="With --hpc: send the run itself to the cluster as a small job (the "
    "profile's 'env' and 'driver') and return; it resumes by itself after a "
    "preemption. Follow it in the run's run.log.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Show detailed logs with timestamps and all model calls. Default: clean progress output.",
)
def run_cmd(
    path: str | None,
    mode: str,
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    n_steps: int | None,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: str,
    hpc: str | None,
    detach: bool,
    verbose: bool,
) -> None:
    """Generate and score candidate molecules across rounds, saving the result.

    The models and folders are the project's; each run writes a new dated folder
    under its rounds/, round by round, so a stopped run can be continued with
    `chemsampler resume`. With --hpc, each round's models run as parallel jobs on
    a cluster (see docs/hpc.md).
    """
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
    # the models it used and how it was launched: enough to run it again
    record = {
        **inputs.describe(),
        "run": {
            "chemsampler_version": _chemsampler_version(),
            "mode": mode,
            "seed_smiles": seed_smiles,
            "original_seed_smiles": original_seed_smiles,
            "n_rounds": n_rounds,
            "n_steps": n_steps,
            "tolerance": tolerance,
            "tanimoto_cutoff": tanimoto_cutoff,
            "tanimoto_direction": tanimoto_direction,
            "hpc": hpc,
            "verbose": verbose,
            "status": "running",
        },
    }
    run_dir, input_dir = inputs.new_round_dir(), Project.load(root).input
    if detach:
        detach_run(inputs, record, run_dir, root, input_dir)
    else:
        execute(inputs, record, run_dir, input_dir)
