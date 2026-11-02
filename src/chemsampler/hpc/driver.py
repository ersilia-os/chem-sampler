"""The run itself as a Slurm job (`run --hpc NAME --detach`, `resume --detach`).

The job only ever runs `chemsampler resume <run>`: its first start, a start after
a preemption (Slurm requeues it) and a resume asked for by hand all continue the
run from the rounds its folder holds. What it prints goes to `run.log` in the run
folder, each start after the last.
"""

import shlex
from pathlib import Path

import click

from .. import store
from .executor import prepare
from .slurm import driver_argv

DRIVER_SCRIPT = "driver.sh"
LOG_FILE = "run.log"

_SCRIPT = """#!/bin/bash
# Written by chemsampler: the run itself, as a job. Slurm starts it again after a
# preemption, and it carries on from the rounds the run folder holds.
set -euo pipefail
unset PYTHONPATH
export PATH={env}/bin:"$PATH"
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export CHEMSAMPLER_HOME={home}
cd {project}
echo "== $(date '+%F %T') job $SLURM_JOB_ID on $(hostname -s), start $((${{SLURM_RESTART_COUNT:-0}} + 1))"
exec chemsampler resume -p {project} {run}
"""


def submit_driver(
    name: str, inputs, run_dir, project_root, *, before_submit=None, say=click.echo
) -> str:
    """
    Check that the run can go to the cluster, then send the run itself as a job.

    Parameters
    ----------
    name : str
        The cluster profile.
    inputs : Inputs
        The run's inputs.
    run_dir : str or Path
        The run's folder.
    project_root : str or Path
        The project the run belongs to.
    before_submit : callable, optional
        Called with the checked `SlurmExecutor` just before the job is sent, e.g.
        to make the run's folder or to check that nothing runs it already.
    say : callable, optional
        Where progress lines go.

    Returns
    -------
    str
        The job id.

    Raises
    ------
    ValueError
        Listing every problem found; nothing has been submitted.
    """
    executor = prepare(name, inputs, run_dir, say=say)
    profile = executor.profile
    wants = profile.driver_resources()
    problems = []
    if profile.env is None:
        problems.append(
            f"profile {name!r} has no 'env': the environment with chemsampler, "
            "visible from the nodes, that the run's own job runs"
        )
    elif not (Path(profile.env) / "bin" / "chemsampler").is_file():
        problems.append(f"no chemsampler in {Path(profile.env) / 'bin'}")
    if not problems:
        # the run's own job may go to other partitions and nodes than the models
        problems += executor.cluster.preflight(
            [Path(profile.env), store.home()],
            partition=wants["partition"],
            nodes=tuple(wants["nodes"] or ()),
        )
    if problems:
        raise ValueError(
            f"cannot send the run to the cluster of profile {name!r}:\n  - "
            + "\n  - ".join(problems)
        )
    if before_submit is not None:
        before_submit(executor)
    q = shlex.quote
    script = Path(run_dir) / DRIVER_SCRIPT
    script.write_text(
        _SCRIPT.format(
            env=q(profile.env),
            home=q(str(store.home().absolute())),
            project=q(str(Path(project_root).absolute())),
            run=q(Path(run_dir).name),
        )
    )
    script.chmod(0o755)
    [job_id] = executor.cluster.submit([driver_argv(profile, run_dir)])
    return job_id


def say_detached(run_dir, job_id: str, say=click.echo) -> None:
    """How to follow, stop and continue a run sent as a job."""
    run_dir = Path(run_dir)
    say(f"Sent the run {run_dir.name} to the cluster as Slurm job {job_id}.")
    say(f"  Follow it:  tail -f {run_dir / LOG_FILE}")
    say(f"  Stop it:    scancel {job_id}")
    say(f"  Continue it once stopped:  chemsampler resume {run_dir.name} --detach")
