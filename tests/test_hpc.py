import json
import signal
import subprocess
import sys
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
    parse_states,
    sbatch_argv,
)
from chemsampler.hub.envs import write_marker
from chemsampler.project import create_project

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
    "source,model,cutoff,direction,priority\n"
    "ersilia,eos4zfy,0,higher,primary\n"
    "rdkit,qed,0-1,range,secondary\n"
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
    # the rdkit environment: its python is the one running these tests
    rdkit = store_root / "envs" / "rdkit"
    (rdkit / "bin").mkdir(parents=True)
    python = rdkit / "bin" / "python"
    python.write_text(f'#!/bin/bash\nexec {sys.executable} "$@"\n')
    python.chmod(0o755)
    write_marker(rdkit, model="rdkit")


@pytest.fixture
def project(tmp_path, store_root, monkeypatch):
    """A project whose three Hub models and RDKit annotator a fake node can run."""
    monkeypatch.setattr(config, "validate_models", lambda **kwargs: None)
    _fake_models(store_root)
    root = tmp_path / "proj"
    create_project(root)
    (root / "input" / "generators.csv").write_text("generators\neos9taz\n")
    (root / "input" / "annotators.csv").write_text(ANNOTATORS)
    profiles = {"test": {**PROFILES["test"], "poll_seconds": 0.01}}
    (store_root / "hpc.json").write_text(json.dumps(profiles))
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
            if self.hold(folder):
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

    def preflight(self, folders):
        self.preflights.append([str(f) for f in folders])
        return list(self.problems)


@pytest.fixture
def cluster(monkeypatch):
    fake = FakeCluster()
    monkeypatch.setattr(executor_module, "Cluster", fake)
    return fake


def _run(root, *args):
    return CliRunner().invoke(
        cli,
        ["run", "-p", str(root), "--mode", "weighted", "--seed-smiles", "CCO"]
        + ["--n-rounds", "2", *args],
    )


def _runs(root):
    return sorted((root / "rounds").iterdir())


# --- a whole run on the "cluster" --------------------------------------------------------


def test_a_run_on_the_cluster_gives_the_same_tables_as_a_local_run(project, cluster):
    local = _run(project)
    on_cluster = _run(project, "--hpc", "test")

    assert local.exit_code == 0, local.output
    assert on_cluster.exit_code == 0, on_cluster.output
    here, there = _runs(project)
    names = sorted(p.name for p in here.glob("*.csv"))
    assert names == ["round0.csv", "round1.csv", "round2.csv", "summary.csv"]
    for name in names:  # every value equal, RDKit's qed included
        pd.testing.assert_frame_equal(
            pd.read_csv(here / name, float_precision="round_trip"),
            pd.read_csv(there / name, float_precision="round_trip"),
            check_exact=True,
        )
    assert json.loads((there / "inputs.json").read_text())["run"]["hpc"] == "test"


def test_each_step_of_a_round_is_one_wave_of_jobs_recorded_in_jobs_json(
    project, cluster
):
    result = _run(project, "--hpc", "test")

    assert result.exit_code == 0, result.output
    [run] = _runs(project)
    record = json.loads((run / "jobs.json").read_text())
    assert record["profile"] == "test"
    waves = {}
    for job in record["jobs"]:
        waves.setdefault(job["stage"], set()).add(job["model"])
        assert job["state"] == "COMPLETED" and job["node"] == "node-a"
        assert (run / job["folder"] / "job.sh").is_file()
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
    assert pd.read_csv(run / "summary.csv")["round"].tolist() == [0, 1]
    record = json.loads((run / "jobs.json").read_text())
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
    assert pd.read_csv(run / "summary.csv")["round"].tolist() == [0, 1]
    record = json.loads((run / "jobs.json").read_text())
    states = {job["job_id"]: job["state"] for job in record["jobs"]}
    assert all(states[i] == "CANCELLED" for i in round2)


def test_problems_found_before_the_run_are_listed_and_nothing_is_sent(project, cluster):
    cluster.problems = ["no node node-b", "the cluster cannot see /somewhere"]

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "cannot run on the cluster of profile 'test'" in result.output
    assert "no node node-b" in result.output
    assert "the cluster cannot see /somewhere" in result.output
    assert cluster.argvs == [] and _runs(project) == []
    [folders] = cluster.preflights  # the rounds, models, envs and rdkit env
    assert any(f.endswith("rounds") for f in folders)
    assert any(f.endswith("rdkit") for f in folders)


def test_a_missing_rdkit_environment_is_found_before_asking_the_cluster(
    project, cluster, store_root
):
    (store_root / "envs" / "rdkit" / ".chemsampler-env.json").unlink()

    result = _run(project, "--hpc", "test")

    assert result.exit_code == 1
    assert "no complete rdkit environment" in result.output
    assert "chemsampler fetch" in result.output
    assert cluster.preflights == [] and cluster.argvs == []


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

    assert profile.ssh == "login" and profile.poll_seconds == 60
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
        ({"poll_seconds": -1}, "'poll_seconds' must be a number"),
        ({"models": {"eos8nyg": {"gpus": 1}}}, r"models.eos8nyg: unknown key"),
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
