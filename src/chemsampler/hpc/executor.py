"""Running a round's models as Slurm jobs, side by side (`run --hpc NAME`)."""

import json
import os
import re
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import click
import pandas as pd

from ..execution import Ran
from ..hub.envs import RDKIT_ENV, env_path, is_complete
from ..models import _rdkit_scores
from ..models.annotator import _RDKitAnnotator
from .jobs import Job, stage_folder, write_model_job, write_rdkit_job
from .profiles import Profile, load_profile
from .slurm import DONE, FAILED, PENDING, RUNNING, Cluster, sbatch_argv

#: How long a job's output may take to show on the shared disk once Slurm says
#: the job completed: network filesystems cache what a folder holds.
OUTPUT_GRACE_SECONDS = 60

#: Jobs are checked every `_FAST_CHECK_SECONDS` for the first `_FAST_FOR_SECONDS`
#: of a step (most jobs take seconds), then every `_SLOW_CHECK_SECONDS`.
_FAST_CHECK_SECONDS = 10
_FAST_FOR_SECONDS = 300
_SLOW_CHECK_SECONDS = 60

#: How long Slurm may not know a submitted job before it counts as lost.
_UNKNOWN_SECONDS = 600

#: How long checks may fail in a row (ssh down, Slurm not answering) before the
#: run stops.
_CHECK_FAILURES_SECONDS = 180

RECORD_FILE = "jobs.json"


def _stage_number(label: str) -> int | None:
    """A step's place in the run, from its label alone: the same on every attempt.

    Round r's sampling is 2r and its annotation 2r + 1 (the seed's is 1).
    """
    found = re.fullmatch(r"Round (\d+) · (sampling|annotation)", label)
    if found is None:
        return None
    return 2 * int(found.group(1)) + (found.group(2) == "annotation")


class SlurmExecutor:
    """
    Send every model call to a Slurm cluster as a job (see `chemsampler.execution`).

    Calls made side by side through `map` become parallel jobs, and each waits
    for its own. All of them are checked together, every 10 s in a step's first
    5 min and every 60 s after, and a line says how many are pending, running and done whenever that changes. A
    job that fails stops the run. Leaving the `with` block on any exception,
    Ctrl-C included, cancels every job still on the cluster. `jobs.json` in the
    run folder records every job.

    Parameters
    ----------
    profile : Profile
        The cluster profile.
    run_dir : str or Path
        The run's folder; the jobs' folders go in its `jobs/`.
    rdkit_env : Path, optional
        The `rdkit` environment, needed when RDKit annotators are used.
    cluster : Cluster, optional
        What submits and checks the jobs (by default, the profile's).
    say : callable, optional
        Where progress lines go.
    """

    def __init__(
        self,
        profile: Profile,
        run_dir,
        *,
        rdkit_env=None,
        cluster=None,
        say=click.echo,
    ):
        self.profile = profile
        self.run_dir = Path(run_dir)
        self.rdkit_env = rdkit_env
        self.cluster = cluster or Cluster(profile)
        self.say = say
        self._lock = threading.Condition()
        self._jobs: list[Job] = []
        self._label, self._number, self._began = "jobs", 0, time.monotonic()
        self._counts = None
        self._failure = None
        self._stop = threading.Event()
        self._watcher = None
        #: What earlier attempts at this run (stopped, then resumed) left in jobs.json.
        self._earlier: list[dict] = []
        self._attempt = 1

    # --- the executor interface (chemsampler.execution) --------------------------

    def map(self, function, items) -> list:
        """`function` applied to each item in its own thread, results in order."""
        items = list(items)
        if len(items) < 2:
            return [function(item) for item in items]
        with ThreadPoolExecutor(max_workers=len(items)) as pool:
            futures = [pool.submit(function, item) for item in items]
            try:
                return [future.result() for future in futures]
            except BaseException:
                self._halt()  # the others stop waiting, so the pool can close
                raise

    def run_model(self, model, smiles_list: list[str]) -> Ran:
        """Run a Hub model as one job and read its output as a local run would."""
        job = self._new_job(model.model_id)
        write_model_job(job, model, smiles_list)
        self._run(job)
        return Ran(pd.read_csv(job.folder / "output.csv"), job.state.elapsed or 0.0)

    def run_rdkit(self, name: str, smiles_list: list[str]) -> Ran:
        """Score with an RDKit annotator as one job, with the `rdkit` environment."""
        if not smiles_list:
            return Ran(pd.DataFrame({"smiles": [], "value": []}, dtype=object), 0.0)
        if self.rdkit_env is None:
            raise RuntimeError(f"{name}: no rdkit environment to run it on the cluster")
        job = self._new_job(name)
        write_rdkit_job(job, name, smiles_list, self.rdkit_env)
        self._run(job)
        values = _rdkit_scores.read(job.folder / "output.csv")
        table = pd.DataFrame({"smiles": smiles_list, "value": values}, dtype=object)
        return Ran(table, job.state.elapsed or 0.0)

    @contextmanager
    def stage(self, label: str):
        """Name a step of a round; its jobs share a folder and a progress line."""
        with self._lock:
            number = _stage_number(label)
            self._number = self._number + 1 if number is None else number
            self._label, self._began, self._counts = label, time.monotonic(), None
            self._set_aside(stage_folder(self.run_dir, self._number, label))
        yield
        jobs = self._stage_jobs()
        if jobs:
            self.say(f"{label}: {len(jobs)} job(s) done ({self._minutes()})")
        self._write_record()

    # --- the run's lifetime ---------------------------------------------------------

    def __enter__(self):
        self._take_over()
        return self

    def job_alive(self, job_id: str) -> bool:
        """Whether a Slurm job is still pending or running (one a run's lock names)."""
        state = self.cluster.states([job_id]).get(job_id)
        return state is not None and state.category in (PENDING, RUNNING)

    def _take_over(self) -> None:
        """Carry on the record of an earlier attempt, cancelling the jobs it left."""
        try:
            data = json.loads((self.run_dir / RECORD_FILE).read_text())
            self._earlier = list(data["jobs"])
        except (OSError, ValueError, KeyError, TypeError):
            return
        self._attempt = 1 + max(r.get("attempt", 1) for r in self._earlier)
        ids = [
            r["job_id"]
            for r in self._earlier
            if r.get("job_id") and r.get("state") != "COMPLETED"
        ]
        if not ids:
            return
        states = self.cluster.states(ids)
        live = [i for i, s in states.items() if s.category in (PENDING, RUNNING)]
        if live:
            self.cluster.cancel(live)
            self.say(f"Cancelled {len(live)} job(s) an earlier attempt left running.")
            states.update(self.cluster.states(live))
        for record in self._earlier:
            if record.get("job_id") in states:
                record["state"] = states[record["job_id"]].state
        self._write_record()

    def _set_aside(self, folder: Path) -> None:
        """Move an earlier attempt's folder for this step out of the way (lock held)."""
        if not folder.exists():
            return
        attempt = f"attempt{self._attempt - 1}"
        target = self.run_dir / "jobs" / "stopped" / attempt / folder.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(folder), str(target))
        old = f"{folder.relative_to(self.run_dir)}{os.sep}"
        new = f"{target.relative_to(self.run_dir)}{os.sep}"
        for record in self._earlier:
            if str(record.get("folder", "")).startswith(old):
                record["folder"] = new + record["folder"][len(old) :]

    def __exit__(self, kind, error, traceback):
        self._halt()
        unfinished = [
            job.job_id
            for job in self._jobs
            if job.job_id
            and (job.state is None or job.state.category in (PENDING, RUNNING))
        ]
        if unfinished:
            live = self._still_live(unfinished)
            ids = unfinished if live is None else live
            try:
                if ids:
                    self.cluster.cancel(ids)
                if live is None:
                    self.say(f"Asked Slurm to cancel {len(ids)} job(s).")
                elif live:
                    self.say(f"Cancelled {len(live)} job(s) still on the cluster.")
                self._last_look(unfinished)
            except RuntimeError as e:
                self.say(f"Could not cancel job(s) {', '.join(ids)}: {e}")
        self._write_record()
        return False

    def _still_live(self, ids: list[str]) -> list[str] | None:
        """
        Those of `ids` Slurm says are pending or running (finished ones are not
        cancelled, nor counted), or None if Slurm cannot be asked.
        """
        try:
            states = self.cluster.states(ids)
        except RuntimeError:
            return None
        # one it does not know may be a job just submitted: cancel it too
        return [
            i
            for i in ids
            if i not in states or states[i].category in (PENDING, RUNNING)
        ]

    def _last_look(self, ids: list[str]) -> None:
        """Ask once more, so `jobs.json` says CANCELLED rather than what was last seen."""
        try:
            states = self.cluster.states(ids)
        except RuntimeError:
            return  # the record keeps what was last known
        for job in self._jobs:
            if job.job_id in states:
                job.state = states[job.job_id]

    # --- jobs -----------------------------------------------------------------------

    def _new_job(self, model_id: str) -> Job:
        with self._lock:
            base = stage_folder(self.run_dir, self._number, self._label)
            folder, n = base / model_id, 2
            while folder.exists():  # a model called twice in one step (chunks)
                folder, n = base / f"{model_id}-{n}", n + 1
            folder.mkdir(parents=True)
            return Job(stage=self._label, model_id=model_id, folder=folder)

    def _run(self, job: Job) -> None:
        argv = sbatch_argv(self.profile, job.model_id, job.name, job.folder)
        [job.job_id] = self.cluster.submit([argv])
        with self._lock:
            self._jobs.append(job)
            self._write_record()
            if self._watcher is None:
                self._watcher = threading.Thread(
                    target=self._watch, name="chemsampler-slurm", daemon=True
                )
                self._watcher.start()
            while True:
                if self._failure is not None:
                    raise RuntimeError(self._failure)
                if self._stop.is_set():
                    raise RuntimeError("the run was stopped")
                if job.state is not None and job.state.category == DONE:
                    break
                self._lock.wait(timeout=1.0)
        self._wait_for_output(job)

    def _wait_for_output(self, job: Job) -> None:
        output = job.folder / "output.csv"
        deadline = time.monotonic() + OUTPUT_GRACE_SECONDS
        while not output.exists():
            if time.monotonic() > deadline:
                raise RuntimeError(
                    f"{job.model_id}: Slurm job {job.job_id} completed but its output "
                    f"never appeared; see {job.folder / 'slurm.err'}"
                )
            if self._stop.wait(1.0):
                raise RuntimeError("the run was stopped")

    def _watch(self) -> None:
        failing_since = None
        while not self._stop.is_set():
            with self._lock:
                live = [
                    j
                    for j in self._jobs
                    if j.state is None or j.state.category in (PENDING, RUNNING)
                ]
            if live:
                try:
                    states = self.cluster.states([j.job_id for j in live])
                    failing_since = None
                except RuntimeError as e:
                    states = None
                    failing_since = failing_since or time.monotonic()
                    if time.monotonic() - failing_since >= _CHECK_FAILURES_SECONDS:
                        self._fail(f"could not check the jobs: {e}")
                with self._lock:
                    if states is not None:
                        self._update(live, states)
                        self._write_record()
                    self._report()
                    self._lock.notify_all()
            self._stop.wait(self._check_interval())

    def _check_interval(self) -> float:
        young = time.monotonic() - self._began < _FAST_FOR_SECONDS
        return _FAST_CHECK_SECONDS if young else _SLOW_CHECK_SECONDS

    def _update(self, live: list[Job], states: dict) -> None:
        """Take in what Slurm said (called with the lock held)."""
        for job in live:
            state = states.get(job.job_id)
            if state is None:
                job.unknown_since = job.unknown_since or time.monotonic()
                if time.monotonic() - job.unknown_since >= _UNKNOWN_SECONDS:
                    self._fail(
                        f"{job.model_id}: Slurm no longer knows job {job.job_id}"
                    )
                continue
            job.unknown_since, job.state = None, state
            if state.category == FAILED:
                self._fail(self._failure_message(job))

    def _fail(self, message: str) -> None:
        with self._lock:
            if self._failure is None:
                self._failure = message
            self._lock.notify_all()

    def _failure_message(self, job: Job) -> str:
        state = job.state
        log = job.folder / "slurm.err"
        on = f" on {state.node}" if state.node else ""
        message = (
            f"{job.model_id}: Slurm job {job.job_id} ended {state.state}{on} "
            f"(exit code {state.exit_code or '?'}); log: {log}"
        )
        try:
            tail = log.read_text(errors="replace").strip().splitlines()[-10:]
        except OSError:
            tail = []
        return "\n".join([message, *tail])

    def _halt(self) -> None:
        self._stop.set()
        with self._lock:
            self._lock.notify_all()

    # --- what the user sees, and what the run folder keeps -------------------------

    def _stage_jobs(self) -> list[Job]:
        return [job for job in self._jobs if job.stage == self._label]

    def _minutes(self) -> str:
        return f"{int((time.monotonic() - self._began) // 60)} min"

    def _report(self) -> None:
        """A progress line when the counts of this step change (lock held)."""
        jobs = self._stage_jobs()
        if not jobs:
            return
        tally = {PENDING: 0, RUNNING: 0, DONE: 0, FAILED: 0}
        for job in jobs:
            tally[PENDING if job.state is None else job.state.category] += 1
        counts = tuple(tally.values())
        if counts == self._counts:
            return
        self._counts = counts
        failed = f", {tally[FAILED]} failed" if tally[FAILED] else ""
        self.say(
            f"{self._label}: {len(jobs)} job(s): {tally[PENDING]} PD, "
            f"{tally[RUNNING]} R, {tally[DONE]} done{failed} ({self._minutes()})"
        )

    def _write_record(self) -> None:
        """`jobs.json`, kept current, so a run killed outright still leaves it."""
        with self._lock:
            jobs = list(self._jobs)
            records = list(self._earlier)
            for job in jobs:
                record = job.record()
                record["folder"] = str(job.folder.relative_to(self.run_dir))
                record["attempt"] = self._attempt
                records.append(record)
        if not records:
            return
        text = json.dumps({"profile": self.profile.name, "jobs": records}, indent=2)
        path = self.run_dir / RECORD_FILE
        path.with_suffix(".partial").write_text(text + "\n")
        os.replace(path.with_suffix(".partial"), path)


def prepare(name: str, inputs, run_dir, *, say=click.echo) -> SlurmExecutor:
    """
    Load profile `name` and check, before any job is sent, that this run can go
    to the cluster: environments for every model (and the `rdkit` one for RDKit
    annotators), Slurm reachable, its partition and nodes there, and the run's
    folders visible from it.

    Raises
    ------
    ValueError
        Listing every problem found; nothing has been submitted.
    """
    profile = load_profile(name)
    problems = []
    envs = None if inputs.envs_dir is None else Path(inputs.envs_dir)
    if envs is None:
        problems.append(
            "the models have no environments folder: on a cluster each model "
            "runs in the environment `chemsampler fetch` builds"
        )
    rdkit_env = None
    if envs is not None and any(
        isinstance(s.annotator, _RDKitAnnotator) for s in inputs.annotators
    ):
        rdkit_env = env_path(envs, RDKIT_ENV)
        if not is_complete(rdkit_env):
            problems.append(
                f"no complete rdkit environment at {rdkit_env}, which RDKit "
                "annotators run with on a cluster; run `chemsampler fetch`"
            )
    cluster = Cluster(profile)
    if not problems:
        folders = [Path(run_dir).parent, inputs.models_dir, envs]
        folders += [] if rdkit_env is None else [rdkit_env]
        try:
            problems += cluster.preflight(folders)
        except RuntimeError as e:
            problems.append(str(e))
    if problems:
        raise ValueError(
            f"cannot run on the cluster of profile {name!r}:\n  - "
            + "\n  - ".join(problems)
        )
    return SlurmExecutor(
        profile, run_dir, rdkit_env=rdkit_env, cluster=cluster, say=say
    )
