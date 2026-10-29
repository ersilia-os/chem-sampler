import os
import subprocess
import tempfile
from typing import Literal

import pandas as pd

from ..utils.logging import logger
from .envs import env_path, is_complete

Backend = Literal["ersilia", "run_sh"]

#: Hub models that look up a seed's nearest neighbours in a reference library
#: instead of generating molecules. Only these can back a `HubSimilarity`, and
#: `HubGenerator` and `HubAnnotator` refuse them.
SIMILARITY_MODEL_IDS = frozenset({"eos8nyg"})


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

    Two backends are available:

    - **"run_sh" (default)**: Shells out directly to the model's bundled run.sh
      script. Requires: model directory at models_dir/model_id with run.sh at
      model/framework/run.sh.
      No persistent server, no port/process contention, much faster startup.
      Avoids Ersilia initialization entirely. Use this unless the model isn't
      conda-packed or isn't fetched locally.

    - **"ersilia"**: Uses Ersilia's serve/run/close API over HTTP. Requires:
      Ersilia installed and the model already fetched via `ersilia fetch`.
      Slower startup (server initialization) but works for any Ersilia model
      (Docker, conda, or other packaging). Use this if the model isn't
      conda-packed or for models fetched via Ersilia directly.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the model (e.g. "eos9taz").
    backend : {"ersilia", "run_sh"}, optional
        Which backend to use, by default "run_sh".
    models_dir : str, optional
        Folder holding one sub-folder per model (`models_dir/model_id`).
        Required by the "run_sh" backend, unused by "ersilia".
    envs_dir : str, optional
        Folder holding one conda environment per model (`envs_dir/model_id`), as
        built by `chemsampler fetch`. With it, the "run_sh" backend runs the
        model inside its own environment, which is what a plain git clone needs.
        Without it the model is expected to bring its own interpreter (a model
        fetched by `ersilia fetch`).

    Raises
    ------
    ValueError
        If `backend` is unknown, or is "run_sh" and `models_dir` is not given.
    """

    def __init__(
        self,
        model_id: str,
        backend: Backend = "run_sh",
        models_dir: str | None = None,
        envs_dir: str | None = None,
    ):
        if backend not in ("ersilia", "run_sh"):
            raise ValueError(f"backend must be 'ersilia' or 'run_sh', got {backend!r}")
        if backend == "run_sh" and not models_dir:
            raise ValueError(
                f"{model_id}: the run_sh backend needs the folder that holds the "
                "model repositories (--path-to-model-repositories)"
            )
        self.model_id = model_id
        self.backend = backend
        self.models_dir = models_dir
        self.envs_dir = envs_dir
        self._ersilia_model = None

    def _get_ersilia_model(self):
        """Lazily import and create ErsiliaModel only when needed."""
        if self._ersilia_model is None:
            from ersilia import ErsiliaModel

            self._ersilia_model = ErsiliaModel(model=self.model_id)
        return self._ersilia_model

    @property
    def model(self):
        """Access the ErsiliaModel instance (lazily loaded on first access)."""
        return self._get_ersilia_model()

    def run(self, smiles_list: list[str]) -> pd.DataFrame:
        """
        Run the model on a list of SMILES and return its raw output table.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to run the model on.

        Returns
        -------
        pandas.DataFrame
            The model's output table, as written to its output CSV.
        """
        if self.backend == "ersilia":
            return self._run_via_ersilia(smiles_list)
        return self._run_via_run_sh(smiles_list)

    def _run_via_ersilia(self, smiles_list: list[str]) -> pd.DataFrame:
        # Ersilia's API is file-only: it rejects an in-memory list, accepts an
        # input CSV of exactly one column, and returns only the path to its
        # output file. Hence the round trip through disk.
        #
        # serve() lives inside the try too: a failed serve() (e.g. a port
        # reused from a previous call) still leaves close() to run, instead of
        # skipping cleanup and propagating the exception straight out of
        # hill_climb with every prior round's result lost.
        model = self._get_ersilia_model()
        try:
            model.serve()
            with tempfile.TemporaryDirectory(prefix="chemsampler-") as tmp_dir:
                input_csv = os.path.join(tmp_dir, "input.csv")
                output_csv = os.path.join(tmp_dir, "output.csv")
                pd.DataFrame({"smiles": smiles_list}).to_csv(input_csv, index=False)
                model.run(input=input_csv, output=output_csv)
                return pd.read_csv(output_csv)
        finally:
            model.close()

    def _process_env(self) -> dict[str, str] | None:
        """
        The environment variables for the model's process, or None to inherit ours.

        A model cloned and built by `chemsampler fetch` has a plain `run.sh`
        (`python main.py ...`), so it must run as if its conda env were active:
        that env's `bin` first on PATH and its own libraries first on the loader
        path (otherwise the system's older libstdc++ is loaded), the user site off,
        and nothing of chemsampler's own PYTHONPATH.
        """
        if not self.envs_dir:
            return None
        env_dir = env_path(os.path.expanduser(self.envs_dir), self.model_id)
        if not is_complete(env_dir):
            raise RuntimeError(
                f"{self.model_id}: no complete environment at {env_dir}; "
                "run `chemsampler fetch` first"
            )
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

    def _run_via_run_sh(self, smiles_list: list[str]) -> pd.DataFrame:
        # No server, no port, no session: run.sh is a one-shot subprocess that
        # bakes its own interpreter path in at pack time, so nothing here needs
        # conda activation, a cwd override, or env changes. Avoids Ersilia entirely.
        bundle = os.path.join(os.path.expanduser(self.models_dir), self.model_id)

        if not os.path.isdir(bundle):
            raise RuntimeError(
                f"{self.model_id}: not found at {bundle}. "
                "Check --path-to-model-repositories, or use backend='ersilia'."
            )

        framework_dir = os.path.join(bundle, "model", "framework")
        run_sh_path = os.path.join(framework_dir, "run.sh")

        if not os.path.isfile(run_sh_path):
            raise RuntimeError(
                f"{self.model_id}: run.sh not found at {run_sh_path}. "
                "Model may not be fully fetched or conda-packed. "
                "Use backend='ersilia' or fetch the model first."
            )
        app_dir = os.path.join(bundle, "app")
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
