import os
import subprocess
import tempfile
from typing import Literal

import pandas as pd
from ersilia import ErsiliaModel

from ..utils.logging import logger

Backend = Literal["ersilia", "run_sh"]


class HubModel:
    """
    Thin wrapper for running an Ersilia Model Hub model on a list of SMILES.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the model (e.g. "eos9taz").
    backend : {"ersilia", "run_sh"}, optional
        How to run the model, by default "ersilia" (serve/run/close over HTTP).
        "run_sh" shells out to the model's bundled `run.sh` directly instead: no
        persistent server, so it sidesteps the orphaned-process and port-contention
        failures the "ersilia" backend can hit under sustained use, at the cost of
        only working for a model that is locally fetched and conda-packed.
    """

    def __init__(self, model_id: str, backend: Backend = "ersilia"):
        if backend not in ("ersilia", "run_sh"):
            raise ValueError(f"backend must be 'ersilia' or 'run_sh', got {backend!r}")
        self.model_id = model_id
        self.backend = backend
        self.model = ErsiliaModel(model=model_id)

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
        try:
            self.model.serve()
            with tempfile.TemporaryDirectory(prefix="chemsampler-") as tmp_dir:
                input_csv = os.path.join(tmp_dir, "input.csv")
                output_csv = os.path.join(tmp_dir, "output.csv")
                pd.DataFrame({"smiles": smiles_list}).to_csv(input_csv, index=False)
                self.model.run(input=input_csv, output=output_csv)
                return pd.read_csv(output_csv)
        finally:
            self.model.close()

    def _run_via_run_sh(self, smiles_list: list[str]) -> pd.DataFrame:
        # No server, no port, no session: run.sh is a one-shot subprocess that
        # bakes its own interpreter path in at pack time, so nothing here needs
        # conda activation, a cwd override, or env changes.
        bundle = self.model.paths["repository"]
        if bundle is None:
            raise RuntimeError(
                f"{self.model_id}: not fetched locally, or fetched but incomplete. "
                "Fetch it first, or use backend='ersilia'."
            )

        service_class_path = os.path.join(bundle, "service_class.txt")
        if not os.path.exists(service_class_path):
            raise RuntimeError(
                f"{self.model_id}: service_class.txt not found under {bundle}. Run it "
                "once with backend='ersilia' so ersilia can determine its packaging, "
                "then retry with backend='run_sh'."
            )
        with open(service_class_path) as f:
            service_class = f.read().strip()
        if service_class != "conda":
            raise RuntimeError(
                f"{self.model_id}: backend='run_sh' only supports conda-packed models "
                f"(found {service_class!r}); use backend='ersilia' instead."
            )

        framework_dir = os.path.join(bundle, "model", "framework")
        run_sh_path = os.path.join(framework_dir, "run.sh")
        app_dir = os.path.join(bundle, "app")

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
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"{self.model_id}: run.sh failed (exit {e.returncode})\n"
                    f"stdout: {e.stdout}\nstderr: {e.stderr}"
                ) from e
            logger.debug(result.stdout)
            return pd.read_csv(output_csv)
