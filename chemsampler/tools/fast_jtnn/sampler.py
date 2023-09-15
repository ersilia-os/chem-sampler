import os
import shutil
import pandas as pd
import tempfile

from ...utils.conda import SimpleConda

root = os.path.dirname(os.path.abspath(__file__))

MOLER_CONDA_ENVIRONMENT = "fastjtnn"


class _JtnnSampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.VOCAB_FILE = os.path.join(root, "data", "vocab.txt")
        self.model_path = os.path.join(
            root, "fast_molvae", "vae_model", "model.epoch-19"
        )
        self.tmp_folder = tempfile.mkdtemp()
        self.data_file = os.path.join(self.tmp_folder, "data_file.txt")
        self.hidden = 450

    def _sample(self, n):
        cmd = "cd {0}; python fast_molvae/sample.py --nsample {1} --vocab {2} --hidden {3} --model {4} --output_file {5}; cd {6}".format(
            self.exec_folder,
            n,
            self.VOCAB_FILE,
            self.hidden,
            self.model_path,
            self.data_file,
            self.cwd,
        )
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)

    def _read_molecules(self):
        with open(self.data_file, "r") as f:
            smiles_list = f.readlines()
        return smiles_list

    def sample(self, n):
        self._sample(n)
        molecules = self._read_molecules()
        return molecules
