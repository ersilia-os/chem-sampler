"""One cluster job: its folder in the run, and the `job.sh` a node runs.

A node needs neither conda nor chemsampler: a Hub model runs its `run.sh` with
its own environment first on PATH (as `HubModel` does here), an RDKit annotator
runs a copy of `_rdkit_scores.py` with the `rdkit` environment's Python. The
output is written under another name and renamed when complete, so a job that is
preempted and starts over never leaves half a file behind.
"""

import re
import shlex
import shutil
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..models import _rdkit_scores
from .slurm import JobState

_HEADER = """#!/bin/bash
# Written by chemsampler: {what}, on {count} molecule(s).
set -euo pipefail
unset PYTHONPATH
export PATH={env}/bin:"$PATH"
export LD_LIBRARY_PATH={env}/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
threads=${{SLURM_CPUS_PER_TASK:-1}}
export OMP_NUM_THREADS=$threads MKL_NUM_THREADS=$threads OPENBLAS_NUM_THREADS=$threads
cd {folder}
rm -f output.csv output.partial.csv
"""


@dataclass
class Job:
    """A model call sent to the cluster, and what Slurm last said about it."""

    stage: str  # e.g. "Round 1 · sampling"
    model_id: str  # a Hub model id, or `mw` / `qed` / `tanimoto`
    folder: Path
    job_id: str | None = None
    state: JobState | None = None
    #: Since when Slurm has not known the job (just submitted, or lost), if it does not.
    unknown_since: float | None = None

    @property
    def name(self) -> str:
        """The Slurm job name, model first (`squeue` shows 8 characters): `eos9p57-r1`."""
        found = re.search(r"Round (\d+)", self.stage)
        return self.model_id if found is None else f"{self.model_id}-r{found.group(1)}"

    def record(self) -> dict:
        """What `jobs.json` keeps about the job."""
        s = self.state
        return {
            "stage": self.stage,
            "model": self.model_id,
            "job_id": self.job_id,
            "state": None if s is None else s.state,
            "node": None if s is None else s.node or None,
            "submitted": None if s is None else s.submitted or None,
            "started": None if s is None else s.started or None,
            "ended": None if s is None else s.ended or None,
            "elapsed_seconds": None if s is None else s.elapsed,
            "exit_code": None if s is None else s.exit_code or None,
            "folder": str(self.folder),
        }


def stage_folder(run_dir, number: int, label: str) -> Path:
    """`jobs/<NN>-<label as a word>`: e.g. `jobs/02-round1-sampling`."""
    word = (
        re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-").replace("round-", "round")
    )
    return Path(run_dir) / "jobs" / f"{number:02d}-{word}"


def write_model_job(job: Job, model, smiles_list: list[str]) -> None:
    """The input and `job.sh` of a Hub model's job (see `HubModel.files`)."""
    framework, run_sh, app = model.files()
    env = model.env_dir()
    if env is None:
        raise RuntimeError(
            f"{model.model_id}: a cluster job needs the environment "
            "`chemsampler fetch` builds; this model has none"
        )
    _write_input(job, smiles_list)
    q = shlex.quote
    script = _HEADER.format(
        what=model.model_id,
        count=len(smiles_list),
        env=q(str(env)),
        folder=q(str(job.folder)),
    )
    script += f"bash {q(run_sh)} {q(framework)} input.csv output.partial.csv {q(app)}\n"
    _finish_script(job, script)


def write_rdkit_job(
    job: Job, name: str, smiles_list: list[str], rdkit_env, args=()
) -> None:
    """The input, scorer and `job.sh` of an in-house RDKit annotator's job."""
    _write_input(job, smiles_list)
    shutil.copy(_rdkit_scores.__file__, job.folder / "_rdkit_scores.py")
    q = shlex.quote
    script = _HEADER.format(
        what=name,
        count=len(smiles_list),
        env=q(str(rdkit_env)),
        folder=q(str(job.folder)),
    )
    extra = "".join(f" {q(arg)}" for arg in args)  # the seed, for tanimoto
    script += f"python _rdkit_scores.py {q(name)} input.csv output.partial.csv{extra}\n"
    _finish_script(job, script)


def _write_input(job: Job, smiles_list: list[str]) -> None:
    job.folder.mkdir(parents=True, exist_ok=True)
    # the same input file a local run writes for run.sh
    pd.DataFrame({"smiles": smiles_list}).to_csv(job.folder / "input.csv", index=False)


def _finish_script(job: Job, script: str) -> None:
    script += "mv output.partial.csv output.csv\n"
    path = job.folder / "job.sh"
    path.write_text(script)
    path.chmod(0o755)
