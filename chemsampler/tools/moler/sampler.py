import os
import shutil
import pandas as pd
import tempfile

from ...utils.conda import SimpleConda

root = os.path.dirname(os.path.abspath(__file__))

MOLER_CONDA_ENVIRONMENT = "moler"

class _MolerSampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.CHECKPOINT_DIR = os.path.join(root, "molecule_generation", "MODEL_DIR" )
        self.tmp_folder = tempfile.mkdtemp()
        self.data_file = os.path.join(self.tmp_folder,"data_file.csv")


    def _sample(self, n):
        cmd = "cd {0}; python ./molecule_generation/cli/sample.py {1} {2} {3} ; cd {4}".format(
            self.exec_folder, self.CHECKPOINT_DIR, n, self.data_file,  self.cwd)
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)
        
    def _read_molecules(self):
        output_file = self.data_file
        smiles_df = pd.read_csv(output_file)
        smiles_list = smiles_df['smiles']
        return smiles_list

    
    def sample(self, n):
        self._sample(n)
        molecules = self._read_molecules()
        return molecules
