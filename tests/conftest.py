import sys

import pytest

from chemsampler import store
from chemsampler.hub.envs import RDKIT_ENV, write_marker
from chemsampler.project import Project, create_project

#: One Hub model of each kind, no RDKit rows and no seed: what tests of other
#: features build on, so they do not change whenever `init`'s templates do.
MINIMAL_INPUTS = {
    "similarity.csv": "model,version\neos8nyg,live\n",
    "generators.csv": "model,version\neos9taz,live\n",
    "annotators.csv": (
        "source,model,version,cutoff,direction,priority\n"
        "ersilia,eos4zfy,live,0.0,higher,primary\n"
    ),
}


def make_project(root, models=None, envs=None):
    """`create_project`, with `MINIMAL_INPUTS` in place of the templates."""
    created = create_project(root, models, envs)
    folder = Project.load(root).input
    (folder / "smiles.csv").unlink()
    for name, text in MINIMAL_INPUTS.items():
        (folder / name).write_text(text)
    return created


def make_rdkit_env(envs_dir):
    """
    A finished `rdkit` environment in `envs_dir`, as `fetch rdkit` leaves it.

    Its `python` is the one running the tests, which has RDKit, so the RDKit
    annotators really run in a process of their own, with this environment.
    """
    env = envs_dir / RDKIT_ENV
    (env / "bin").mkdir(parents=True, exist_ok=True)
    python = env / "bin" / "python"
    python.write_text(f'#!/bin/bash\nexec {sys.executable} "$@"\n')
    python.chmod(0o755)
    write_marker(env, model=RDKIT_ENV)
    return env


@pytest.fixture(autouse=True)
def store_root(tmp_path, monkeypatch):
    """
    Point chemsampler's folder somewhere harmless for every test.

    Without this a test would read, and `fetch` would write, the real
    `~/.chemsampler` (or `$CHEMSAMPLER_HOME`) of whoever runs the suite.
    """
    root = tmp_path / "store"
    monkeypatch.setenv(store.HOME_VARIABLE, str(root))
    # a suite run inside a Slurm job must not look like a cluster run
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    return root


@pytest.fixture
def fake_fetch(store_root):
    """
    Make models look fetched into the test store, as `chemsampler fetch` leaves
    them: a clone folder and an environment, finished unless `complete` is False.
    """

    def fetch(*model_ids, complete=True):
        for model_id in model_ids:
            (store_root / "models" / model_id).mkdir(parents=True)
            env = store_root / "envs" / model_id
            env.mkdir(parents=True)
            if complete:
                write_marker(env, model=model_id)

    return fetch
