import os
import tempfile

import pandas as pd
from ersilia import ErsiliaModel


class HubModel:
    """
    Thin wrapper for running an Ersilia Model Hub model on a list of SMILES.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the model (e.g. "eos9taz").
    """

    def __init__(self, model_id: str):
        self.model_id = model_id
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
        self.model.serve()
        tmp_dir = tempfile.mkdtemp(prefix="chemsampler-")
        input_csv = os.path.join(tmp_dir, "input.csv")
        output_csv = os.path.join(tmp_dir, "output.csv")
        pd.DataFrame({"smiles": smiles_list}).to_csv(input_csv, index=False)
        self.model.run(input=input_csv, output=output_csv)
        df = pd.read_csv(output_csv)
        self.model.close()
        return df
