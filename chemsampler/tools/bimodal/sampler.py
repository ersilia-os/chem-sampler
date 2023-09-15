from ...utils.conda import SimpleConda
import os
import pandas as pd
import tempfile


root = os.path.dirname(os.path.abspath(__file__))

BIMODAL_CONDA_ENVIRONMENT = "bimodal"


class _BimodalSampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.tmp_folder = tempfile.mkdtemp()
        self.data_file = os.path.join(self.tmp_folder, "data_file.csv")

    def _sample(self, n):
        cmd = "cd {0} ; python ./model/_sampler.py {1} {2}; cd {3}".format(
            self.exec_folder, n, self.data_file, self.cwd
        )
        SimpleConda().run_commandlines(BIMODAL_CONDA_ENVIRONMENT, cmd)

    def _read_molecules(self):
        output_file = self.data_file
        smiles_df = pd.read_csv(output_file)
        smiles_list = smiles_df["smiles"]
        return smiles_list

    def sample(self, n):
        print(self.cwd)
        print(self.exec_folder)
        self._sample(n)
        molecules = self._read_molecules()
        return molecules
