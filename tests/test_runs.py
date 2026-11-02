import json
import socket
import subprocess

import pytest

from chemsampler import runs


def _held_by(run_dir, **who):
    run_dir.mkdir(exist_ok=True)
    held = {"host": socket.gethostname(), "pid": None, "job": None, **who}
    (run_dir / runs.LOCK_FILE).write_text(json.dumps(held))


def _dead_pid():
    done = subprocess.Popen(["true"])
    done.wait()
    return done.pid


def test_a_free_run_is_taken(tmp_path):
    runs.take_lock(tmp_path)

    held = json.loads((tmp_path / runs.LOCK_FILE).read_text())
    assert held["host"] == socket.gethostname()
    runs.release_lock(tmp_path)
    assert not (tmp_path / runs.LOCK_FILE).exists()


def test_a_lock_left_by_a_dead_process_here_is_taken_over(tmp_path):
    _held_by(tmp_path, pid=_dead_pid())

    runs.take_lock(tmp_path)


def test_a_run_held_on_another_machine_is_refused(tmp_path):
    _held_by(tmp_path, host="elsewhere", pid=1)

    with pytest.raises(ValueError, match="may still be going on elsewhere.*delete"):
        runs.take_lock(tmp_path)


def test_a_run_held_by_a_live_job_is_refused_and_a_finished_one_taken(tmp_path):
    _held_by(tmp_path, host="node-a", job="42")

    with pytest.raises(ValueError, match="still going as Slurm job 42.*scancel 42"):
        runs.take_lock(tmp_path, job_alive=lambda job: True)
    runs.take_lock(tmp_path, job_alive=lambda job: False)


def test_a_job_cannot_be_checked_without_its_cluster(tmp_path):
    _held_by(tmp_path, host="node-a", job="42")

    with pytest.raises(ValueError, match="may still be going as Slurm job 42"):
        runs.take_lock(tmp_path)


def test_a_requeued_job_takes_its_own_run_back(tmp_path, monkeypatch):
    # preempted, then started again by Slurm: same job id, any node
    monkeypatch.setenv("SLURM_JOB_ID", "42")
    _held_by(tmp_path, host="node-a", job="42")

    runs.take_lock(tmp_path, job_alive=lambda job: True)


def test_the_latest_run_or_a_named_one_is_found(tmp_path):
    for name in ("2026-11-01_10-00-00", "2026-11-02_09-00-00", "notes"):
        (tmp_path / name).mkdir()
    for name in ("2026-11-01_10-00-00", "2026-11-02_09-00-00"):
        (tmp_path / name / runs.RECORD_FILE).write_text("{}")

    assert runs.find_run(tmp_path).name == "2026-11-02_09-00-00"
    assert runs.find_run(tmp_path, "2026-11-01_10-00-00").name == (
        "2026-11-01_10-00-00"
    )
    with pytest.raises(ValueError, match="notes is not a run"):
        runs.find_run(tmp_path, "notes")
