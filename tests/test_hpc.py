import json
import signal
import subprocess
import threading
import time
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from chemsampler import config
from chemsampler.cli import cli
from chemsampler.hpc import executor as executor_module
from chemsampler.hpc.jobs import Job
from chemsampler.hpc.profiles import Profile, load_profile
from chemsampler.hpc.slurm import (
    DONE,
    FAILED,
    PENDING,
    RUNNING,
    Cluster,
    JobState,
    driver_argv,
    parse_states,
    sbatch_argv,
)
from chemsampler.hub.envs import write_marker

from .conftest import make_project, make_rdkit_env

# --- fake models that a "node" (bash) can run -----------------------------------------

#: run.sh bodies: generator and similarity write fixed neighbours; the annotator
#: scores each SMILES by its length, so values are easy to predict.
RUN_SH = {
    "eos9taz": "printf 'smi_0,smi_1\\nCCN,CCC\\n' > \"$3\"\n",
    "eos8nyg": "printf 'smi_00\\nCCCl\\n' > \"$3\"\n",
    "eos4zfy": (
        'awk -F, \'NR==1{print "smiles,score"} NR>1{print $1","length($1)}\' '
        '"$2" > "$3"\n'
    ),
}
ANNOTATORS = (
    "source,model,version,cutoff,direction,priority\n"
    "ersilia,eos4zfy,v1.0.0,1000,higher,primary\n"
    "rdkit,qed,,0-1,range,guard-rail\n"
)
PROFILES = {"test": {"partition": "fake", "nodes": ["node-a", "node-b"]}}


def _fake_models(store_root, run_sh=RUN_SH):
    for model, body in run_sh.items():
        framework = store_root / "models" / model / "model" / "framework"
        framework.mkdir(parents=True)
        (framework / "run.sh").write_text(body)
        env = store_root / "envs" / model
        (env / "bin").mkdir(parents=True)
        write_marker(env, model=model)
    make_rdkit_env(store_root / "envs")


@pytest.fixture
def project(tmp_path, store_root, monkeypatch):
    """A project whose three Hub models and RDKit annotator a fake node can run."""
    monkeypatch.setattr(config, "validate_models", lambda **kwargs: None)
    _fake_models(store_root)
    root = tmp_path / "proj"
    make_project(root)
    (root / "input" / "generators.csv").write_text("model,version\neos9taz,v1.0.0\n")
    (root / "input" / "annotators.csv").write_text(ANNOTATORS)
    (store_root / "hpc.json").write_text(json.dumps(PROFILES))
    monkeypatch.setattr(executor_module, "_FAST_CHECK_SECONDS", 0.01)
    monkeypatch.setattr(executor_module, "_SLOW_CHECK_SECONDS", 0.01)
    return root


class FakeCluster:
    """
    Stands in for Slurm: a job runs (bash job.sh, as a node would) the first time
    it is checked. `hold(folder)` keeps a job pending instead; `on_check` sees
    every check first.
    """

    def __init__(self):
        self.argvs, self.folders, self.outcomes = [], {}, {}
        self.cancelled, self.preflights, self.problems = [], [], []
        self.hold = lambda folder: False
        self.queued = set()  # job ids that stay pending until cancelled
        self.placements = []
        self.on_check = None

    def __call__(self, profile):  # used in place of the Cluster class
        self.profile = profile
        return self

    def submit(self, argvs):
        ids = []
        for argv in argvs:
            job_id = str(1000 + len(self.argvs))
            self.argvs.append(argv)
            self.folders[job_id] = Path(argv[-1]).parent
            ids.append(job_id)
        return ids

    def states(self, ids):
        if self.on_check is not None:
            self.on_check(self, ids)
        found = {}
        for job_id in ids:
            folder = self.folders[job_id]
            if job_id in self.cancelled:  # sacct keeps it, as Slurm does
                found[job_id] = JobState(job_id, "CANCELLED", FAILED)
                continue
            if job_id in self.queued or self.hold(folder):
                found[job_id] = JobState(job_id, "PENDING", PENDING)
                continue
            if job_id not in self.outcomes:
                ran = subprocess.run(
                    ["bash", str(folder / "job.sh")],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                (folder / "slurm.err").write_text(ran.stderr)
                ok = ran.returncode == 0
                self.outcomes[job_id] = JobState(
                    job_id,
                    "COMPLETED" if ok else "FAILED",
                    DONE if ok else FAILED,
                    node="node-a",
                    elapsed=2.0,
                    exit_code=f"{ran.returncode}:0",
                )
            found[job_id] = self.outcomes[job_id]
        return found

    def cancel(self, ids):
        self.cancelled.extend(ids)

    def preflight(self, folders, partition=None, nodes=None):
        self.preflights.append([str(f) for f in folders])
        self.placements.append((partition, nodes))
        return list(self.problems)


@pytest.fixture
def cluster(monkeypatch):
    fake = FakeCluster()
    monkeypatch.setattr(executor_module, "Cluster", fake)
    return fake


def _run(root, *args):
    return CliRunner().invoke(
        cli,
        ["run", "-p", str(root), "--seed-smiles", "CCO"] + ["--n-rounds", "2", *args],
    )


def _runs(root):
    return sorted((root / "runs").iterdir())


# --- a whole run on the "cluster" --------------------------------------------------------


def test_a_run_on_the_cluster_gives_the_same_tables_as_a_local_run(project, cluster):
    local = _run(project)
    on_cluster = _run(project, "--hpc", "test")

    assert local.exit_code == 0, local.output
    assert on_cluster.exit_code == 0, on_cluster.output
    here, there = _runs(project)
    names = sorted(p.name for p in (here / "seed1").glob("*.csv"))
    assert names == ["round0.csv", "round1.csv", "round2.csv", "summary.csv"]
    for name in names:  # every value equal, RDKit's qed included
        pd.testing.assert_frame_equal(
            pd.read_csv(here / "seed1" / name, float_precision="round_trip"),
            pd.read_csv(there / "seed1" / name, float_precision="round_trip"),
            check_exact=True,
        )
    assert json.loads((there / "inputs.json").read_text())["run"]["hpc"] == "test"


def test_each_step_of_a_round_is_one_wave_of_jobs_recorded_in_jobs_json(
    project, cluster
):
    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    [run] = _runs(project)
    record = json.loads((run / "seed1" / "jobs.json").read_text())
    assert record["profile"] == "test"
    waves = {}
    for job in record["jobs"]:
        waves.setdefault(job["stage"], set()).add(job["model"])
        assert job["state"] == "COMPLETED" and job["node"] == "node-a"
        assert (run / "seed1" / job["folder"] / "job.sh").is_file()
    assert waves == {
        "Round 0 · annotation": {"eos4zfy", "qed"},
        "Round 1 · sampling": {"eos9taz", "eos8nyg"},
        "Round 1 · annotation": {"eos4zfy", "qed"},
        "Round 2 · sampling": {"eos9taz", "eos8nyg"},
        "Round 2 · annotation": {"eos4zfy", "qed"},
    }
    assert "Round 1 · sampling: 2 job(s) done" in result.output
    for argv in cluster.argvs:  # only where the profile allows
        assert "--partition=fake" in argv and "--nodelist=node-a,node-b" in argv


def test_each_seed_s_jobs_and_jobs_json_are_in_its_own_folder(project, cluster):
    (project / "input" / "smiles.csv").write_text("smiles\nCCO\nCCCO\n")

    result = CliRunner().invoke(
        cli, ["run", "-p", str(project), "--n-rounds", "1", "--hpc", "test"]
    )

    assert result.exit_code == 0, result.output
    [run] = _runs(project)
    for seed in ("seed1", "seed2"):
        record = json.loads((run / seed / "jobs.json").read_text())
        assert len(record["jobs"]) == 6  # 2 annotators, then 2 samplers + 2
        assert all((run / seed / job["folder"]).is_dir() for job in record["jobs"])
    assert not (run / "jobs").exists() and not (run / "jobs.json").exists()
    assert len(cluster.argvs) == 12


def test_the_models_of_a_step_wait_on_the_cluster_together(project, cluster):
    # a job is held until its whole step is submitted: run one model at a time and
    # this never happens, the checks run out and the run fails
    checks = []

    def whole_step_submitted(fake, ids):
        checks.append(ids)
        if len(checks) > 2000:
            raise RuntimeError("the step was never submitted as a whole")

    def held(folder):
        siblings = [f for f in cluster.folders.values() if f.parent == folder.parent]
        return len(siblings) < 2

    cluster.on_check, cluster.hold = whole_step_submitted, held

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output


def test_a_failed_job_stops_the_run_cancels_the_rest_and_keeps_finished_rounds(
    project, cluster, store_root
):
    fail_twice = (
        'n=$(cat "$1/count" 2>/dev/null || echo 0); echo $((n+1)) > "$1/count"\n'
        '[ "$n" -ge 1 ] && { echo "boom: out of memory" >&2; exit 1; }\n'
    )
    run_sh = store_root / "models" / "eos9taz" / "model" / "framework" / "run.sh"
    run_sh.write_text(fail_twice + RUN_SH["eos9taz"])
    cluster.hold = lambda f: "round2" in str(f) and f.name == "eos8nyg"

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "eos9taz: Slurm job" in result.output and "FAILED" in result.output
    assert "slurm.err" in result.output and "boom: out of memory" in result.output
    [held] = [i for i, f in cluster.folders.items() if cluster.hold(f)]
    assert cluster.cancelled == [held]
    assert "Cancelled 1 job(s) still on the cluster." in result.output
    [run] = _runs(project)
    assert pd.read_csv(run / "seed1" / "summary.csv")["round"].tolist() == [0, 1]
    record = json.loads((run / "seed1" / "jobs.json").read_text())
    states = {job["job_id"]: job["state"] for job in record["jobs"]}
    assert states[held] == "CANCELLED"


def test_ctrl_c_cancels_the_jobs_on_the_cluster_and_keeps_finished_rounds(
    project, cluster
):
    cluster.hold = lambda f: "round2" in str(f)
    pressed = []

    def interrupt_in_round_2(fake, ids):
        # a real SIGINT, as a terminal sends: it also wakes a thread that waits
        if not pressed and any(fake.hold(fake.folders[i]) for i in ids):
            pressed.append(True)
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    cluster.on_check = interrupt_in_round_2

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    assert "Partial results written" in result.output
    round2 = sorted(i for i, f in cluster.folders.items() if "round2" in str(f))
    assert round2 and sorted(cluster.cancelled) == round2
    [run] = _runs(project)
    assert pd.read_csv(run / "seed1" / "summary.csv")["round"].tolist() == [0, 1]
    record = json.loads((run / "seed1" / "jobs.json").read_text())
    states = {job["job_id"]: job["state"] for job in record["jobs"]}
    assert all(states[i] == "CANCELLED" for i in round2)


def test_ctrl_c_cancels_and_counts_only_the_jobs_slurm_still_runs(
    project, cluster, monkeypatch
):
    pressed, released = [], []
    cluster.hold = lambda f: (
        "round2" in str(f) and not (released and f.name == "eos9taz")
    )

    def interrupt_then_finish_one(fake, ids):
        if pressed:  # the next check: the eos9taz job has finished meanwhile
            released.append(True)
        elif any(fake.hold(fake.folders[i]) for i in ids):
            # the watcher waits long now, so the next check is the one made
            # before cancelling
            monkeypatch.setattr(executor_module, "_FAST_CHECK_SECONDS", 60)
            monkeypatch.setattr(executor_module, "_SLOW_CHECK_SECONDS", 60)
            pressed.append(True)
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    cluster.on_check = interrupt_then_finish_one

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    assert "Cancelled 1 job(s) still on the cluster." in result.output
    round2 = {f.name: i for i, f in cluster.folders.items() if "round2" in str(f)}
    assert cluster.cancelled == [round2["eos8nyg"]]
    [run] = _runs(project)
    record = json.loads((run / "seed1" / "jobs.json").read_text())
    states = {job["job_id"]: job["state"] for job in record["jobs"]}
    assert states[round2["eos9taz"]] == "COMPLETED"
    assert states[round2["eos8nyg"]] == "CANCELLED"


def test_ctrl_c_after_every_job_finished_cancels_nothing(project, cluster, monkeypatch):
    pressed, released = [], []
    cluster.hold = lambda f: "round2" in str(f) and not released

    def interrupt_then_finish_all(fake, ids):
        if pressed:  # the check before cancelling: both jobs finished meanwhile
            released.append(True)
        elif any(fake.hold(fake.folders[i]) for i in ids):
            monkeypatch.setattr(executor_module, "_FAST_CHECK_SECONDS", 60)
            monkeypatch.setattr(executor_module, "_SLOW_CHECK_SECONDS", 60)
            pressed.append(True)
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    cluster.on_check = interrupt_then_finish_all

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    assert cluster.cancelled == []
    assert "Cancelled" not in result.output and "Asked Slurm" not in result.output


def test_when_slurm_cannot_say_what_still_runs_every_unfinished_job_is_cancelled(
    project, cluster, monkeypatch
):
    cluster.hold = lambda f: "round2" in str(f)
    pressed = []

    def interrupt_then_fail(fake, ids):
        if pressed:
            raise RuntimeError("ssh: connection closed")
        if any(fake.hold(fake.folders[i]) for i in ids):
            monkeypatch.setattr(executor_module, "_FAST_CHECK_SECONDS", 60)
            monkeypatch.setattr(executor_module, "_SLOW_CHECK_SECONDS", 60)
            pressed.append(True)
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    cluster.on_check = interrupt_then_fail

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    assert "Asked Slurm to cancel 2 job(s)." in result.output
    round2 = sorted(i for i, f in cluster.folders.items() if "round2" in str(f))
    assert sorted(cluster.cancelled) == round2


def test_checks_that_keep_failing_stop_the_run_after_a_while(
    project, cluster, monkeypatch
):
    monkeypatch.setattr(executor_module, "_CHECK_FAILURES_SECONDS", 0.2)

    def down(fake, ids):
        raise RuntimeError("ssh: connection refused")

    cluster.on_check = down

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "could not check the jobs: ssh: connection refused" in result.output


def test_jobs_are_checked_often_early_in_a_step_then_less(tmp_path):
    executor = executor_module.SlurmExecutor(
        Profile("t"), tmp_path, cluster=FakeCluster()
    )
    assert executor._check_interval() == 10
    executor._began -= 299
    assert executor._check_interval() == 10
    executor._began -= 2
    assert executor._check_interval() == 60


def test_a_job_slurm_does_not_know_counts_as_lost_only_after_ten_minutes(tmp_path):
    executor = executor_module.SlurmExecutor(
        Profile("t"), tmp_path, cluster=FakeCluster()
    )
    job = Job(stage="Round 1 · sampling", model_id="eos9taz", folder=tmp_path)
    job.job_id = "1000"
    executor._update([job], {})
    assert executor._failure is None and job.unknown_since is not None
    job.unknown_since -= 599
    executor._update([job], {})
    assert executor._failure is None
    job.unknown_since -= 2
    executor._update([job], {})
    assert "Slurm no longer knows job 1000" in executor._failure


def test_a_resumed_cluster_run_cancels_what_a_killed_attempt_left_and_keeps_its_record(
    project, cluster
):
    whole = _run(project, "--hpc", "test")
    assert whole.exit_code == 0, whole.output
    # an attempt stopped in round 2 ...
    cluster.hold = lambda f: "round2" in str(f)
    cluster.on_check = lambda fake, ids: (
        any(fake.hold(fake.folders[i]) for i in ids)
        and signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)
    )
    assert _run(project, "--hpc", "test").exit_code == 0
    _, run = _runs(project)
    round2 = sorted(
        i
        for i, f in cluster.folders.items()
        if "round2" in str(f) and f.is_relative_to(run)
    )
    # ... then killed outright: its round 2 jobs still queued, its lock left
    cluster.cancelled.clear()
    record = json.loads((run / "seed1" / "jobs.json").read_text())
    for job in record["jobs"]:
        if job["job_id"] in round2:
            job["state"] = "PENDING"
    (run / "seed1" / "jobs.json").write_text(json.dumps(record))
    cluster.hold, cluster.on_check = (lambda f: False), None
    cluster.queued = set(round2)

    result = CliRunner().invoke(cli, ["resume", "-p", str(project), run.name])

    assert result.exit_code == 0, result.output
    assert "Cancelled 2 job(s) an earlier attempt left running." in result.output
    assert sorted(cluster.cancelled) == round2
    first, _ = _runs(project)
    for name in ("round0.csv", "round1.csv", "round2.csv", "summary.csv"):
        pd.testing.assert_frame_equal(
            pd.read_csv(first / "seed1" / name),
            pd.read_csv(run / "seed1" / name),
            check_exact=True,
        )
    jobs = json.loads((run / "seed1" / "jobs.json").read_text())["jobs"]
    assert {job["attempt"] for job in jobs} == {1, 2}
    redone = [job for job in jobs if job["attempt"] == 2]
    assert {job["stage"] for job in redone} == {
        "Round 2 · sampling",
        "Round 2 · annotation",
    }  # rounds 0 and 1 were read back, not sent again
    stale = run / "seed1" / "jobs" / "stopped" / "attempt1" / "04-round2-sampling"
    assert stale.is_dir() and (run / "seed1" / "jobs" / "04-round2-sampling").is_dir()
    assert all(
        job["folder"].startswith("jobs/stopped/attempt1/")
        for job in jobs
        if job["job_id"] in round2
    )


def test_jobs_json_lists_a_job_while_it_waits(project, cluster):
    # a run killed outright (preempted) still leaves the ids of its jobs
    listed = []

    def look(fake, ids):
        path = _runs(project)[0] / "seed1" / "jobs.json"
        record = json.loads(path.read_text()) if path.exists() else {"jobs": []}
        listed.append(set(ids) <= {job["job_id"] for job in record["jobs"]})

    cluster.on_check = look

    assert _run(project, "--hpc", "test").exit_code == 0
    assert listed and all(listed)


# --- the run itself as a job (--detach) ------------------------------------------------


def _with_env(project, store_root, **driver):
    """Give the test profile an environment with chemsampler (and driver settings)."""
    env = store_root / "chemsampler-env"
    (env / "bin").mkdir(parents=True)
    (env / "bin" / "chemsampler").write_text("#!/bin/bash\n")
    profiles = json.loads((store_root / "hpc.json").read_text())
    profiles["test"].update(env=str(env), driver=driver)
    (store_root / "hpc.json").write_text(json.dumps(profiles))
    return env


def _resume(project, *args, **env):
    return CliRunner().invoke(cli, ["resume", "-p", str(project), *args], env=env)


def test_a_detached_run_is_sent_as_one_job_that_resumes_it(
    project, cluster, store_root
):
    env = _with_env(project, store_root, partition="other", nodes=["node-c"])

    result = _run(project, "--hpc", "test", "--detach")

    assert result.exit_code == 0, result.output
    [run] = _runs(project)
    [argv] = cluster.argvs  # only the run's own job: no model ran here
    assert argv[-1] == f"{run}/driver.sh" and "--requeue" in argv
    [job_id] = cluster.folders
    assert f"Slurm job {job_id}" in result.output
    assert f"tail -f {run / 'run.log'}" in result.output
    assert f"scancel {job_id}" in result.output
    assert f"chemsampler resume {run.name} --detach" in result.output
    record = json.loads((run / "inputs.json").read_text())["run"]
    assert record["status"] == "queued" and record["driver_job"] == job_id
    assert (run / "input" / "annotators.csv").is_file()
    assert not list(run.glob("round*.csv"))
    script = (run / "driver.sh").read_text()
    assert f'export PATH={env}/bin:"$PATH"' in script
    assert f"export CHEMSAMPLER_HOME={store_root}" in script
    assert f"exec chemsampler resume -p {project} {run.name}" in script
    # the cluster sees the env and the store, on the partition the job goes to
    assert [str(env), str(store_root)] in cluster.preflights
    assert ("other", ("node-c",)) in cluster.placements


def test_the_run_s_own_job_runs_it_where_slurm_is(project, cluster, store_root):
    _with_env(project, store_root)
    assert _run(project, "--hpc", "test", "--detach").exit_code == 0
    [run] = _runs(project)
    [job_id] = cluster.folders
    cluster.queued.add(job_id)  # it is running: it is the job below

    result = _resume(project, run.name, SLURM_JOB_ID=job_id)

    assert result.exit_code == 0, result.output
    assert "Resuming" in result.output and "(queued)" in result.output
    assert json.loads((run / "inputs.json").read_text())["run"]["status"] == "finished"
    assert (run / "seed1" / "summary.csv").is_file()


def test_a_run_waiting_as_a_job_is_not_run_twice(project, cluster, store_root):
    _with_env(project, store_root)
    assert _run(project, "--hpc", "test", "--detach").exit_code == 0
    [run] = _runs(project)
    [job_id] = cluster.folders
    cluster.queued.add(job_id)  # still in the queue: no lock yet

    here = _resume(project, run.name)
    again = _resume(project, run.name, "--detach")

    for result in (here, again):
        assert result.exit_code == 1
        assert f"queued or running as Slurm job {job_id}" in result.output
    assert len(cluster.argvs) == 1
    cluster.cancelled.append(job_id)  # scancel: now it can be sent again
    assert _resume(project, run.name, "--detach").exit_code == 0
    assert len(cluster.argvs) == 2


def test_a_new_run_whose_job_is_refused_leaves_no_folder(
    project, cluster, store_root, monkeypatch
):
    _with_env(project, store_root)

    def refuse(argvs):
        raise RuntimeError("sbatch: error: Time limit specification required")

    monkeypatch.setattr(cluster, "submit", refuse)

    result = _run(project, "--hpc", "test", "--detach")

    assert result.exit_code == 1 and "Time limit specification" in result.output
    assert not list((project / "runs").iterdir())


def test_detaching_needs_a_cluster_and_an_env(project, cluster, store_root):
    without_hpc = _run(project, "--detach")
    without_env = _run(project, "--hpc", "test", "--detach")

    assert without_hpc.exit_code == 1 and "it needs --hpc NAME" in without_hpc.output
    assert without_env.exit_code == 1 and "has no 'env'" in without_env.output
    assert cluster.argvs == [] and not list((project / "runs").iterdir())


def test_problems_found_before_the_run_are_listed_and_nothing_is_sent(project, cluster):
    cluster.problems = ["no node node-b", "the cluster cannot see /somewhere"]

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "cannot run on the cluster of profile 'test'" in result.output
    assert "no node node-b" in result.output
    assert "the cluster cannot see /somewhere" in result.output
    assert cluster.argvs == [] and _runs(project) == []
    [folders] = cluster.preflights  # the runs, models, envs and rdkit env
    assert any(f.endswith("runs") for f in folders)
    assert any(f.endswith("rdkit") for f in folders)


def test_a_missing_rdkit_environment_is_found_before_asking_the_cluster(
    project, cluster, store_root
):
    (store_root / "envs" / "rdkit" / ".chemsampler-env.json").unlink()

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "qed run in the rdkit environment, which is not built" in result.output
    assert "`chemsampler fetch rdkit`" in result.output
    assert cluster.preflights == [] and cluster.argvs == []


def test_a_tanimoto_annotator_is_a_job_given_the_seed_with_local_values(
    project, cluster, store_root
):
    (project / "input" / "annotators.csv").write_text(
        "source,model,version,cutoff,direction,priority\n"
        "ersilia,eos4zfy,v1.0.0,1000,higher,primary\n"
        "rdkit,tanimoto,,0-1,range,guard-rail\n"
    )

    local = _run(project)
    result = _run(project, "--hpc", "test")

    assert local.exit_code == 0 and result.exit_code == 0, result.output
    folders = [Path(argv[-1]).parent for argv in cluster.argvs]
    tanimoto = [f for f in folders if f.name == "tanimoto"]
    assert tanimoto and all(
        (f / "job.sh")
        .read_text()
        .rstrip()
        .endswith("output.partial.csv CCO\nmv output.partial.csv output.csv")
        for f in tanimoto
    )
    [preflight] = cluster.preflights
    assert str(store_root / "envs" / "rdkit") in preflight
    here, there = _runs(project)
    for name in ("round1.csv", "round2.csv"):
        pd.testing.assert_frame_equal(
            pd.read_csv(here / "seed1" / name, float_precision="round_trip"),
            pd.read_csv(there / "seed1" / name, float_precision="round_trip"),
            check_exact=True,
        )


def test_an_unknown_profile_says_which_ones_exist(project, cluster):
    result = _run(project, "--hpc", "nope")

    assert result.exit_code == 1
    assert "has no profile 'nope'; it has: test" in result.output


# --- profiles ----------------------------------------------------------------------------


def _write_profiles(store_root, profiles):
    store_root.mkdir(parents=True, exist_ok=True)
    (store_root / "hpc.json").write_text(json.dumps(profiles))


def test_a_profile_gives_each_model_its_own_resources(store_root):
    _write_profiles(
        store_root,
        {
            "c": {
                "ssh": "login",
                "cpus": 4,
                "mem": "8G",
                "time": "01:00:00",
                "models": {"eos8nyg": {"mem": "16G"}},
            }
        },
    )

    profile = load_profile("c")

    assert profile.ssh == "login"
    assert profile.resources("eos8nyg") == {"cpus": 4, "mem": "16G", "time": "01:00:00"}
    assert profile.resources("qed") == {"cpus": 4, "mem": "8G", "time": "01:00:00"}


def test_without_a_profiles_file_it_says_where_to_write_one(store_root):
    with pytest.raises(ValueError, match=f"write {store_root / 'hpc.json'}"):
        load_profile("c")


@pytest.mark.parametrize(
    "profile, message",
    [
        ({"partitions": "x"}, r"unknown key\(s\) \['partitions'\]"),
        ({"nodes": "node-a"}, "'nodes' must be a list of node names"),
        ({"cpus": 0}, "'cpus' must be a whole number"),
        ({"mem": 8}, "'mem' must be text"),
        ({"poll_seconds": 10}, r"unknown key\(s\) \['poll_seconds'\]"),
        ({"models": {"eos8nyg": {"gpus": 1}}}, r"models.eos8nyg: unknown key"),
        ({"env": 3}, "'env' must be text"),
        ({"driver": {"gpus": 1}}, r"driver: unknown key"),
        ({"driver": {"cpus": 0}}, "driver: 'cpus' must be a whole number"),
        ({"driver": {"nodes": "node-a"}}, "driver: 'nodes' must be a list"),
    ],
)
def test_a_malformed_profile_is_explained(store_root, profile, message):
    _write_profiles(store_root, {"c": profile})

    with pytest.raises(ValueError, match=message):
        load_profile("c")


# --- what is sent to Slurm, and what is read back -----------------------------------------


def test_sbatch_asks_for_the_profile_s_nodes_and_the_model_s_resources():
    profile = Profile(
        "c",
        partition="spot",
        qos="low",
        nodes=("node-a", "node-b"),
        cpus=2,
        mem="4G",
        models={"eos8nyg": {"cpus": 8, "time": "02:00:00"}},
    )

    argv = sbatch_argv(profile, "eos8nyg", "eos8nyg-r1", "/runs/x")

    assert argv[:2] == ["sbatch", "--parsable"] and argv[-1] == "/runs/x/job.sh"
    for expected in [
        "--job-name=eos8nyg-r1",
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=8",
        "--requeue",
        "--partition=spot",
        "--qos=low",
        "--nodelist=node-a,node-b",
        "--mem=4G",
        "--time=02:00:00",
        "--output=/runs/x/slurm.out",
        "--error=/runs/x/slurm.err",
    ]:
        assert expected in argv
    plain = sbatch_argv(Profile("p"), "qed", "qed-r1", "/runs/y")
    assert not [a for a in plain if a.split("=")[0] in ("--nodelist", "--partition")]


def test_slurm_states_are_read_from_squeue_then_sacct():
    text = (
        "11|PD|\n"
        "12|R|node-a\n"
        "--\n"
        "11|PENDING|None assigned|0|0:0|2026-11-01T10:00:00|Unknown|Unknown\n"
        "12|RUNNING|node-a|30|0:0|2026-11-01T10:00:00|2026-11-01T10:01:00|Unknown\n"
        "13|COMPLETED|node-b|95|0:0|2026-11-01T10:00:00|2026-11-01T10:01:00|"
        "2026-11-01T10:02:35\n"
        "13.batch|COMPLETED|node-b|95|0:0|x|x|x\n"
        "14|OUT_OF_MEMORY|node-a|12|0:125|x|x|x\n"
        "15|CANCELLED by 1234|node-a|5|0:0|x|x|x\n"
        "16|PREEMPTED|node-a|5|0:0|x|x|x\n"  # left the queue: it is not coming back
    )

    states = parse_states(text)

    assert {i: s.category for i, s in states.items()} == {
        "11": PENDING,
        "12": RUNNING,
        "13": DONE,
        "14": FAILED,
        "15": FAILED,
        "16": FAILED,
    }
    assert states["11"].node == "" and states["13"].elapsed == 95.0
    assert (
        states["13"].ended == "2026-11-01T10:02:35"
        and states["15"].state == "CANCELLED"
    )


def test_a_requeued_job_back_in_the_queue_counts_as_pending():
    states = parse_states("17|PD|\n--\n17|PREEMPTED|node-a|40|0:0|x|x|x\n")

    assert states["17"].category == PENDING


def test_commands_go_over_ssh_when_the_profile_names_a_host(monkeypatch):
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs["input"]))
        return subprocess.CompletedProcess(argv, 0, "4242;cluster\n", "")

    monkeypatch.setattr("chemsampler.hpc.slurm.subprocess.run", run)

    ids = Cluster(Profile("c", ssh="login")).submit([["sbatch", "--parsable", "j.sh"]])
    Cluster(Profile("here")).cancel(["4242"])

    assert ids == ["4242"]
    (over_ssh, script), (here, cancel) = calls
    assert over_ssh[0] == "ssh" and "BatchMode=yes" in over_ssh
    assert over_ssh[-2:] == ["login", "bash -s"] and "sbatch --parsable j.sh" in script
    assert here == ["bash", "-s"] and cancel.strip() == "scancel 4242"


def test_the_checks_can_look_at_another_partition_and_nodes(monkeypatch):
    scripts = []

    def run(argv, **kwargs):
        scripts.append(kwargs["input"])
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr("chemsampler.hpc.slurm.subprocess.run", run)
    cluster = Cluster(Profile("c", partition="spot", nodes=("node-a",)))

    cluster.preflight(["/x"])
    cluster.preflight(["/x"], partition="other", nodes=("node-c",))

    profiles, driver = scripts
    assert "-p spot" in profiles and "-n node-a" in profiles
    assert "-p other" in driver and "-n node-c" in driver
    assert "spot" not in driver and "node-a" not in driver


def test_inside_a_cluster_job_slurm_is_never_reached_over_ssh(monkeypatch):
    # the run's own job (--detach) uses the same profile as the machine it came from
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr("chemsampler.hpc.slurm.subprocess.run", run)
    monkeypatch.setenv("SLURM_JOB_ID", "77")

    Cluster(Profile("c", ssh="login")).cancel(["4242"])

    assert calls == [["bash", "-s"]]


def test_the_run_s_own_job_asks_for_its_resources_and_restarts_after_preemption():
    profile = Profile(
        "c",
        partition="spot",
        qos="low",
        nodes=("node-a",),
        cpus=4,
        mem="8G",
        time="01:00:00",
        driver={"mem": "4G", "nodes": ("node-b",)},
    )

    argv = driver_argv(profile, "/runs/x")

    assert argv[:2] == ["sbatch", "--parsable"] and argv[-1] == "/runs/x/driver.sh"
    for expected in [
        "--job-name=chemsampler",
        "--cpus-per-task=1",  # not the models' 4
        "--requeue",
        "--open-mode=append",
        "--output=/runs/x/run.log",
        "--error=/runs/x/run.log",
        "--partition=spot",
        "--qos=low",
        "--mem=4G",
        "--nodelist=node-b",
        "--time=01:00:00",  # the profile's, unless driver.time says otherwise
    ]:
        assert expected in argv
    longer = Profile("c", time="01:00:00", driver={"time": "3-00:00:00"})
    assert "--time=3-00:00:00" in driver_argv(longer, "/runs/x")


def test_a_job_that_stops_halfway_leaves_no_output(tmp_path):
    from chemsampler.hpc.jobs import write_model_job
    from chemsampler.hub.client import HubModel

    _fake_models(tmp_path, {"eos9taz": 'echo smi_0 > "$3"\nexit 1\n'})
    model = HubModel(
        "eos9taz", models_dir=str(tmp_path / "models"), envs_dir=str(tmp_path / "envs")
    )
    job = Job("Round 1 · sampling", "eos9taz", tmp_path / "job")
    write_model_job(job, model, ["CCO"])

    ran = subprocess.run(["bash", str(job.folder / "job.sh")], check=False)

    # a preempted or failed attempt never leaves a file that looks finished
    assert ran.returncode != 0 and not (job.folder / "output.csv").exists()


# --- a completed job's output on a network disk --------------------------------------------


def _waiting(tmp_path):
    executor = executor_module.SlurmExecutor(
        Profile("c"), tmp_path, cluster=FakeCluster()
    )
    return executor, Job("Round 1 · sampling", "eos9taz", tmp_path, job_id="1")


def test_an_output_that_shows_late_is_waited_for(tmp_path):
    executor, job = _waiting(tmp_path)
    late = threading.Timer(0.3, (tmp_path / "output.csv").write_text, ["smi_0\nCCN\n"])
    late.start()

    executor._wait_for_output(job)

    assert (tmp_path / "output.csv").exists()


def test_an_output_that_never_shows_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr(executor_module, "OUTPUT_GRACE_SECONDS", 0.2)
    executor, job = _waiting(tmp_path)
    start = time.monotonic()

    with pytest.raises(RuntimeError, match="completed but its output never appeared"):
        executor._wait_for_output(job)

    assert time.monotonic() - start < 5
