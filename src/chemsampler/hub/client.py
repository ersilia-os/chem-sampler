import os
import subprocess
import tempfile

import pandas as pd

from .. import execution
from ..utils.logging import logger
from .envs import env_path, is_complete

#: Hub models that look up a seed's nearest neighbours in a reference library
#: instead of generating molecules. Only these can back a `HubSimilarity`, and
#: `HubGenerator` and `HubAnnotator` refuse them. All three run the same search
#: (100 neighbours, ranked by fusing CLAMP, SAND and Mol-JEPA) on a fixed local
#: library: the Ersilia Reference Library, COCONUT natural products and the
#: Enamine Hit Locator. The Hub's other similarity models call outside services.
SIMILARITY_MODEL_IDS = frozenset({"eos8nyg", "eos1bgt", "eos80al"})


def not_allowed_similarity_message(model_id: str) -> str:
    """Error text for a similarity model chemsampler does not support."""
    allowed = ", ".join(sorted(SIMILARITY_MODEL_IDS))
    return (
        f"{model_id} is not meant to be used in chemsampler as a similarity "
        f"model. Choose from: {allowed}."
    )


class HubModel:
    """
    Thin wrapper for running an Ersilia Model Hub model on a list of SMILES.

    Runs the model's own `model/framework/run.sh` from its folder as a one-shot
    subprocess: no server, no port, no Ersilia session.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the model (e.g. "eos9taz").
    models_dir : str
        Folder holding one sub-folder per model (`models_dir/model_id`).
    envs_dir : str, optional
        Folder holding one conda environment per model (`envs_dir/model_id`), as
        built by `chemsampler fetch`. With it, the model runs inside its own
        environment, which is what a plain git clone needs. Without it the model
        is expected to bring its own interpreter (a model fetched by
        `ersilia fetch`).

    Raises
    ------
    ValueError
        If `models_dir` is not given.
    """

    def __init__(
        self,
        model_id: str,
        models_dir: str | None,
        envs_dir: str | None = None,
    ):
        if not models_dir:
            raise ValueError(
                f"{model_id}: needs the folder that holds the model repositories "
                "(models_dir)"
            )
        self.model_id = model_id
        self.models_dir = models_dir
        self.envs_dir = envs_dir
        #: How long the model itself ran on its last call, in seconds.
        self.elapsed = None

    def run(self, smiles_list: list[str]) -> pd.DataFrame:
        """
        Run the model on a list of SMILES and return its raw output table.

        The call goes through the current executor (`chemsampler.execution`):
        here by default, or as a cluster job.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to run the model on.

        Returns
        -------
        pandas.DataFrame
            The model's output table, as written to its output CSV.
        """
        ran = execution.current().run_model(self, smiles_list)
        self.elapsed = ran.elapsed
        return ran.table

    def env_dir(self):
        """
        The model's own conda environment, as `chemsampler fetch` built it, or
        None for a model that brings its own interpreter (`ersilia fetch`).

        Raises
        ------
        RuntimeError
            If an environments folder is set but this model's is not complete.
        """
        if not self.envs_dir:
            return None
        env_dir = env_path(os.path.expanduser(self.envs_dir), self.model_id)
        if not is_complete(env_dir):
            raise RuntimeError(
                f"{self.model_id}: no complete environment at {env_dir}; "
                "run `chemsampler fetch` first"
            )
        return env_dir

    def _process_env(self) -> dict[str, str] | None:
        """
        The environment variables for the model's process, or None to inherit ours.

        A model cloned and built by `chemsampler fetch` has a plain `run.sh`
        (`python main.py ...`), so it must run as if its conda env were active:
        that env's `bin` first on PATH and its own libraries first on the loader
        path (otherwise the system's older libstdc++ is loaded), the user site off,
        and nothing of chemsampler's own PYTHONPATH. A cluster job sets the same
        (`chemsampler.hpc.jobs`).
        """
        env_dir = self.env_dir()
        if env_dir is None:
            return None
        env = {**os.environ}
        env.pop("PYTHONPATH", None)
        env["PATH"] = os.pathsep.join(
            filter(None, [str(env_dir / "bin"), env.get("PATH")])
        )
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            filter(None, [str(env_dir / "lib"), env.get("LD_LIBRARY_PATH")])
        )
        env["PYTHONNOUSERSITE"] = "1"
        env["PYTHONUNBUFFERED"] = "1"
        return env

    def files(self) -> tuple[str, str, str]:
        """
        The model's framework folder, its `run.sh` and its app folder: what
        `bash run.sh <framework> <input.csv> <output.csv> <app>` needs.

        Raises
        ------
        RuntimeError
            If the model or its `run.sh` is not there.
        """
        bundle = os.path.join(os.path.expanduser(self.models_dir), self.model_id)
        if not os.path.isdir(bundle):
            raise RuntimeError(
                f"{self.model_id}: not found at {bundle}; run `chemsampler fetch`"
            )
        framework_dir = os.path.join(bundle, "model", "framework")
        run_sh_path = os.path.join(framework_dir, "run.sh")
        if not os.path.isfile(run_sh_path):
            raise RuntimeError(
                f"{self.model_id}: run.sh not found at {run_sh_path}; the model "
                "is not fully fetched: run `chemsampler fetch`"
            )
        return framework_dir, run_sh_path, os.path.join(bundle, "app")

    def run_here(self, smiles_list: list[str]) -> pd.DataFrame:
        """Run the model on this machine (what `LocalExecutor` does)."""
        # No server, no port, no session: run.sh is a one-shot subprocess, run in
        # the model's own environment (see _process_env) or, for a model fetched
        # by `ersilia fetch`, with the interpreter it baked in at pack time.
        framework_dir, run_sh_path, app_dir = self.files()
        process_env = self._process_env()

        with tempfile.TemporaryDirectory(prefix="chemsampler-") as tmp_dir:
            input_csv = os.path.join(tmp_dir, "input.csv")
            output_csv = os.path.join(tmp_dir, "output.csv")
            pd.DataFrame({"smiles": smiles_list}).to_csv(input_csv, index=False)
            try:
                result = subprocess.run(
                    [
                        "bash",
                        run_sh_path,
                        framework_dir,
                        input_csv,
                        output_csv,
                        app_dir,
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    env=process_env,
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"{self.model_id}: run.sh failed (exit {e.returncode})\n"
                    f"stdout: {e.stdout}\nstderr: {e.stderr}"
                ) from e
            logger.debug(result.stdout)
            return pd.read_csv(output_csv)
