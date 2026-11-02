"""Talking to Slurm: `sbatch`, `squeue`, `sacct` and `scancel`, here or over ssh.

Each call is one small bash script, run with `bash -s` here, or on the login
node a profile names with `ssh <host> bash -s`. Reading what Slurm answers is
done by pure functions, so it can be tested without a cluster.
"""

import os
import shlex
import subprocess
from dataclasses import dataclass

from .profiles import Profile

PENDING, RUNNING, DONE, FAILED = "pending", "running", "done", "failed"

#: `squeue` short states of a job still on its way to finishing.
_RUNNING_CODES = {"R", "CG"}

#: `sacct` states (first word) of a job that has not ended; anything not here or
#: COMPLETED means it ended without its result.
_SACCT_LIVE = {
    "PENDING": PENDING,
    "REQUEUED": PENDING,
    "RUNNING": RUNNING,
    "COMPLETING": RUNNING,
}

_SACCT_FIELDS = "JobID,State,NodeList,ElapsedRaw,ExitCode,Submit,Start,End"


@dataclass(frozen=True)
class JobState:
    """What Slurm says about one job."""

    job_id: str
    state: str  # Slurm's own word: PENDING, RUNNING, COMPLETED, TIMEOUT, ...
    category: str  # pending, running, done or failed
    node: str = ""
    elapsed: float | None = None  # seconds the job has run
    exit_code: str = ""
    submitted: str = ""
    started: str = ""
    ended: str = ""


def sbatch_argv(profile: Profile, model_id: str, name: str, folder) -> list[str]:
    """The `sbatch` command for one job whose `job.sh` is in `folder`."""
    wants = profile.resources(model_id)
    argv = [
        "sbatch",
        "--parsable",
        f"--job-name={name}",
        "--nodes=1",
        "--ntasks=1",
        f"--cpus-per-task={wants['cpus']}",
        # a preempted job goes back to the queue; job.sh starts over cleanly
        "--requeue",
        f"--output={folder}/slurm.out",
        f"--error={folder}/slurm.err",
    ]
    where = {"partition": profile.partition, "qos": profile.qos, "nodes": profile.nodes}
    argv += _placement({**where, **wants})
    argv.append(f"{folder}/job.sh")
    return argv


def driver_argv(profile: Profile, run_dir) -> list[str]:
    """The `sbatch` command for the run's own job (`driver.sh` in `run_dir`)."""
    log = f"{run_dir}/run.log"
    argv = [
        "sbatch",
        "--parsable",
        "--job-name=chemsampler",
        "--nodes=1",
        "--ntasks=1",
        f"--cpus-per-task={profile.driver_resources()['cpus']}",
        # preempted, it goes back to the queue and resumes the run where it was
        "--requeue",
        # every start writes after the last one, not over it
        "--open-mode=append",
        f"--output={log}",
        f"--error={log}",
    ]
    argv += _placement(profile.driver_resources())
    argv.append(f"{run_dir}/driver.sh")
    return argv


def _placement(wants: dict) -> list[str]:
    argv = []
    for key in ("partition", "qos", "mem", "time"):
        if wants[key]:
            argv.append(f"--{key}={wants[key]}")
    if wants["nodes"]:
        argv.append(f"--nodelist={','.join(wants['nodes'])}")
    return argv


def parse_states(text: str) -> dict[str, JobState]:
    """
    Read what `Cluster.states` printed: `squeue` lines (`id|state|nodes`), a
    line `--`, then `sacct -P` lines. A job still in the queue is judged by
    `squeue` (a requeued job is back there), one that has left it by `sacct`.
    """
    lines = [line.strip() for line in text.splitlines()]
    cut = lines.index("--") if "--" in lines else len(lines)
    queued = {}
    for line in lines[:cut]:
        parts = line.split("|")
        if len(parts) == 3:
            queued[parts[0]] = (parts[1], parts[2])
    accounted = {}
    for line in lines[cut + 1 :]:
        parts = line.split("|")
        if len(parts) == 8 and "." not in parts[0]:  # allocations, not steps
            accounted[parts[0]] = parts
    states = {}
    for job_id in {*queued, *accounted}:
        _, word, node, elapsed, exit_code, submitted, started, ended = accounted.get(
            job_id, [job_id, "", "", "", "", "", "", ""]
        )
        word = word.split(" ")[0]  # "CANCELLED by 1234"
        if job_id in queued:
            code, queued_node = queued[job_id]
            category = RUNNING if code in _RUNNING_CODES else PENDING
            if word not in _SACCT_LIVE:  # sacct lags behind, or saw an old attempt
                word = code
            node = queued_node or node
        elif word == "COMPLETED":
            category = DONE
        else:
            category = _SACCT_LIVE.get(word, FAILED)
        states[job_id] = JobState(
            job_id=job_id,
            state=word,
            category=category,
            node="" if node in ("None assigned", "(null)") else node,
            elapsed=float(elapsed) if elapsed.isdigit() else None,
            exit_code=exit_code,
            submitted=_when(submitted),
            started=_when(started),
            ended=_when(ended),
        )
    return states


def _when(value: str) -> str:
    return "" if value in ("Unknown", "None", "") else value


class Cluster:
    """The Slurm commands of one profile, run here or on its login node."""

    def __init__(self, profile: Profile):
        self.profile = profile

    @property
    def here(self) -> bool:
        """Whether Slurm is on this machine: no `ssh`, or inside a cluster job."""
        return self.profile.ssh is None or bool(os.environ.get("SLURM_JOB_ID"))

    @property
    def where(self) -> str:
        return "this machine" if self.here else self.profile.ssh

    def shell(self, script: str) -> str:
        """Run a bash script where Slurm is; return what it printed."""
        if self.here:
            argv = ["bash", "-s"]
        else:
            argv = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30"]
            argv += [self.profile.ssh, "bash -s"]
        try:
            done = subprocess.run(
                argv, input=script, capture_output=True, text=True, check=False
            )
        except OSError as e:
            raise RuntimeError(f"could not reach Slurm on {self.where}: {e}") from e
        if done.returncode != 0:
            said = (done.stderr or done.stdout).strip()[-400:]
            raise RuntimeError(f"Slurm on {self.where} failed: {said}")
        return done.stdout

    def submit(self, argvs: list[list[str]]) -> list[str]:
        """Submit jobs (their `sbatch_argv`); return their ids, in order."""
        out = self.shell("set -e\n" + "\n".join(shlex.join(a) for a in argvs) + "\n")
        ids = [line.split(";")[0].strip() for line in out.splitlines() if line.strip()]
        if len(ids) != len(argvs):
            raise RuntimeError(f"sbatch on {self.where} answered {out.strip()!r}")
        return ids

    def states(self, job_ids: list[str]) -> dict[str, JobState]:
        """What Slurm says about these jobs; a job it does not know is left out."""
        ids = shlex.quote(",".join(job_ids))
        script = (
            f"squeue -h -j {ids} -o '%i|%t|%N' 2>/dev/null || true\n"
            "echo --\n"
            f"sacct -n -X -P -j {ids} -o {_SACCT_FIELDS} 2>/dev/null || true\n"
        )
        return parse_states(self.shell(script))

    def cancel(self, job_ids: list[str]) -> None:
        self.shell("scancel " + " ".join(shlex.quote(i) for i in job_ids) + "\n")

    def preflight(self, folders, partition=None, nodes=None) -> list[str]:
        """
        Problems that would stop jobs from running, one line each (none: ready):
        missing Slurm commands, an unknown partition or node (the profile's,
        unless others are given), a folder the cluster cannot see.
        """
        partition = partition or self.profile.partition
        nodes = self.profile.nodes if nodes is None else nodes
        lines = [
            "for c in sbatch squeue sacct scancel sinfo; do",
            '  command -v "$c" >/dev/null || echo "Slurm command not found: $c"',
            "done",
        ]
        if partition:
            p = shlex.quote(partition)
            lines.append(
                f'[ -n "$(sinfo -h -p {p} -o %P 2>/dev/null)" ] '
                f"|| echo 'no partition {partition}'"
            )
        for node in nodes:
            n = shlex.quote(node)
            lines.append(
                f"[ -n \"$(sinfo -h -n {n} -o %N 2>/dev/null)\" ] || echo 'no node {node}'"
            )
        for folder in folders:
            f = shlex.quote(str(folder))
            lines.append(
                f"[ -d {f} ] || echo {shlex.quote(f'the cluster cannot see {folder}')}"
            )
        return [
            line for line in self.shell("\n".join(lines) + "\n").splitlines() if line
        ]
