"""A run's folder: what the run started with, how far it got, and who runs it.

`chemsampler run` makes the folder before the first round, with a copy of the
inputs and an `inputs.json` of the options. Every round's table is saved in it
as soon as it is made (see `hill_climb`'s `checkpoint`), so a run stopped in any
way, even killed, continues with `chemsampler resume` from what is there.
"""

import json
import os
import shutil
import socket
from pathlib import Path

from .utils.jsonfile import read_object

RECORD_FILE = "inputs.json"

#: The inputs as the run found them; `resume` reads these, not the project's,
#: which may have changed since.
INPUT_FOLDER = "input"
INPUT_FILES = ("similarity.csv", "generators.csv", "annotators.csv", "smiles.csv")

#: Who is running the run right now: host, process and Slurm job.
LOCK_FILE = "run.lock"


def start(run_dir, input_dir, record: dict) -> None:
    """Make a new run's folder: a copy of its inputs, and its record."""
    copy = Path(run_dir) / INPUT_FOLDER
    copy.mkdir(parents=True)
    for name in INPUT_FILES:
        if (Path(input_dir) / name).is_file():
            shutil.copy(Path(input_dir) / name, copy / name)
    write_record(run_dir, record)


def read_record(run_dir) -> dict:
    """The run's `inputs.json`: its models and options, and its `status`."""
    return read_object(Path(run_dir) / RECORD_FILE)


def write_record(run_dir, record: dict) -> None:
    path = Path(run_dir) / RECORD_FILE
    path.with_suffix(".partial").write_text(json.dumps(record, indent=2) + "\n")
    os.replace(path.with_suffix(".partial"), path)


def set_status(run_dir, record: dict, status: str) -> None:
    """Record how the run stands: running, finished, stopped or failed."""
    record["run"]["status"] = status
    write_record(run_dir, record)


def find_run(runs_dir, name: str | None = None) -> Path:
    """
    A run's folder: `name` as a path or as a folder of `runs_dir`, or the latest run.

    Raises
    ------
    ValueError
        If there is no such run (or no run at all).
    """
    folder = Path(runs_dir)
    if name is None:
        runs = []
        if folder.is_dir():
            runs = sorted(p for p in folder.iterdir() if (p / RECORD_FILE).is_file())
        if not runs:
            raise ValueError(f"no runs in {folder} yet")
        return runs[-1]  # the folders are named by the time they started
    for path in (Path(name), folder / name):
        if (path / RECORD_FILE).is_file():
            return Path(os.path.abspath(path))
    raise ValueError(f"{name} is not a run (no {RECORD_FILE} in it, nor in {folder})")


def take_lock(run_dir, job_alive=None) -> None:
    """
    Mark the run as being run here, unless it may still be running elsewhere.

    Parameters
    ----------
    run_dir : str or Path
        The run's folder.
    job_alive : callable, optional
        Whether a Slurm job id is still pending or running; without it, a run
        held by a job is taken to be possibly alive.

    Raises
    ------
    ValueError
        If the lock names a live job (other than this one), a live process on
        this machine, or a process on another machine.
    """
    path = Path(run_dir) / LOCK_FILE
    here = _whoami()
    if path.exists():
        check_lock(run_dir, job_alive)
    path.write_text(json.dumps(here) + "\n")


def check_lock(run_dir, job_alive=None) -> None:
    """Raise ValueError if the run may be going elsewhere (see `take_lock`)."""
    path = Path(run_dir) / LOCK_FILE
    try:
        held = json.loads(path.read_text())
    except FileNotFoundError:
        return
    except (OSError, ValueError):
        held = {}
    here = _whoami()
    job, host, pid = held.get("job"), held.get("host"), held.get("pid")
    if job:
        if job == here["job"]:
            return  # a requeued job, taking its own run back
        if job_alive is None:
            raise ValueError(
                f"the run may still be going as Slurm job {job}; if it is not, "
                f"delete {path}"
            )
        if job_alive(job):
            raise ValueError(
                f"the run is still going as Slurm job {job}; stop it with "
                f"`scancel {job}` first"
            )
        return
    if host == here["host"]:
        if _alive(pid):
            raise ValueError(
                f"the run is still going on this machine (process {pid}); stop it "
                f"first, or delete {path} if that process is not chemsampler"
            )
        return
    raise ValueError(
        f"the run may still be going on {host} (process {pid}); if it is not, "
        f"delete {path}"
    )


def check_driver(record: dict, job_alive) -> None:
    """
    Refuse to run a run that its own cluster job (`--detach`) still runs or waits to.

    A job still in the queue holds no lock yet, so the lock alone cannot tell.

    Raises
    ------
    ValueError
        If `record` names a driver job, other than this one, that is pending or
        running.
    """
    job = record["run"].get("driver_job")
    if job and job != os.environ.get("SLURM_JOB_ID") and job_alive(job):
        raise ValueError(
            f"the run is queued or running as Slurm job {job}; stop it with "
            f"`scancel {job}` first"
        )


def release_lock(run_dir) -> None:
    (Path(run_dir) / LOCK_FILE).unlink(missing_ok=True)


def _whoami() -> dict:
    return {
        "host": socket.gethostname(),
        "pid": os.getpid(),
        "job": os.environ.get("SLURM_JOB_ID"),
    }


def _alive(pid) -> bool:
    if not isinstance(pid, int):
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:  # someone else's process
        return True
    return True
