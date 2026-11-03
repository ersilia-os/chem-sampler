"""Where models run: every model call goes through the current executor.

By default (`LocalExecutor`) each call runs here, one after another, as it is
made. `chemsampler.hpc.SlurmExecutor` instead sends every call to a cluster as a
job: the models of one step of a round (`stage`) are called in parallel through
`map`, so their jobs run at the same time.
"""

import contextlib
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class Ran:
    """What a model call returned, and how long the model itself ran."""

    table: pd.DataFrame
    elapsed: float


class LocalExecutor:
    """Run every model here, one call at a time (chemsampler's default)."""

    def map(self, function, items) -> list:
        """`function` applied to each item, in order, one after another."""
        return [function(item) for item in items]

    def run_model(self, model, smiles_list: list[str]) -> Ran:
        """Run a Hub model's `run.sh` on the SMILES (see `HubModel`)."""
        start = time.monotonic()
        table = model.run_here(smiles_list)
        return Ran(table, time.monotonic() - start)

    def run_rdkit(self, annotator, smiles_list: list[str]) -> Ran:
        """
        Score the SMILES with an in-house RDKit annotator (column `value`): a
        copy of `_rdkit_scores.py` run by the `rdkit` environment's Python, as a
        cluster job runs it, or, without an environments folder, in this process.
        """
        from .hub.envs import process_env
        from .models import _rdkit_scores

        args, env_dir = annotator.args(), annotator.env_dir()
        start = time.monotonic()
        if env_dir is None:
            values = _rdkit_scores.score(annotator.name, smiles_list, *args)
        else:
            with tempfile.TemporaryDirectory(prefix="chemsampler-") as tmp:
                folder = Path(tmp)
                # its own folder: nothing beside the script can shadow an import
                shutil.copy(_rdkit_scores.__file__, folder / "_rdkit_scores.py")
                pd.DataFrame({"smiles": smiles_list}).to_csv(
                    folder / "input.csv", index=False
                )
                command = [str(env_dir / "bin" / "python"), "_rdkit_scores.py"]
                command += [annotator.name, "input.csv", "output.csv", *args]
                ran = subprocess.run(
                    command,
                    cwd=folder,
                    env=process_env(env_dir),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if ran.returncode != 0:
                    raise RuntimeError(
                        f"{annotator.name}: the rdkit environment failed (exit "
                        f"{ran.returncode})\nstderr: {ran.stderr}"
                    )
                values = _rdkit_scores.read(folder / "output.csv")
        table = pd.DataFrame({"smiles": smiles_list, "value": values}, dtype=object)
        return Ran(table, time.monotonic() - start)

    @contextlib.contextmanager
    def stage(self, label: str):
        """Name a step of a round; nothing to do when models run here."""
        yield


_LOCAL = LocalExecutor()
# A plain module setting, not a context variable: the worker threads of
# `SlurmExecutor.map` must see the executor that started them.
_current = None


def current():
    """The executor model calls go through now."""
    return _current or _LOCAL


@contextlib.contextmanager
def use(executor):
    """Send every model call made inside the block through `executor`."""
    global _current
    previous, _current = _current, executor
    try:
        yield executor
    finally:
        _current = previous


def stage(label: str):
    """Name the step of a round the model calls inside the block belong to."""
    return current().stage(label)
