import os
import shutil
import configparser
import tempfile
import pandas as pd
import numpy as np

from ...utils.conda import SimpleConda

root = os.path.dirname(os.path.abspath(__file__))

MOLER_CONDA_ENVIRONMENT = "moler"


class MolerFit:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.CHECKPOINT_PATH = os.path.join(
            root,
            "molecule_generation",
            "MODEL_DIR",
            "GNN_Edge_MLP_MoLeR__2022-02-24_07-16-23_best.pkl",
        )
        self.tuned_checkpoint_dir = os.path.join(root, "molecule_generation", "output")
        self.INPUT_DIR = tempfile.mkdtemp()
        self.OUTPUT_DIR = tempfile.mkdtemp()
        self.TRACE_DIR = tempfile.mkdtemp()
        self.tmp_folder = tempfile.mkdtemp()
        self.data_file = os.path.join(self.tmp_folder, "data_file.csv")
        self.train_split = 0.8
        self.test_split = 0.1
        self.val_split = 0.1

    def input_split(self, smiles_list):
        assert (self.train_split + self.test_split + self.val_split) == 1
        assert self.val_split == self.test_split
        df = pd.DataFrame(smiles_list)

        indices = np.arange(df.shape[0])
        np.random.shuffle(indices)
        train_size = int(df.shape[0] * self.train_split)
        test_size = int(df.shape[0] * self.test_split)
        val_size = df.shape[0] - train_size - test_size

        # split the dataframe using the indices and sizes
        train_df, test_df, val_df = np.split(
            df.iloc[indices], [train_size, train_size + test_size]
        )

        # indices_or_sections = [int(self.train_split * len(df)), int((1 - self.val_split - self.test_split) * len(df))]
        # train_ds, val_ds, test_ds = np.split(df_sample, indices_or_sections)
        # print(len(train_ds), len(val_ds), len(test_ds))
        # return train_ds, val_ds, test_ds
        return train_df, val_df, test_df

    def write_smi(self, smiles, filename):
        csv_filename = filename + ".csv"
        csv_filename = os.path.join(self.INPUT_DIR, csv_filename)
        smi_filename = filename + ".smiles"
        smi_filename = os.path.join(self.INPUT_DIR, smi_filename)
        smiles.to_csv(csv_filename, header=False, index=False)
        os.rename(csv_filename, smi_filename)

    def save_data(self, train_ds, val_ds, test_ds):
        self.write_smi(train_ds, "train")
        self.write_smi(val_ds, "valid")
        self.write_smi(test_ds, "test")

    def _fit(self, smiles_list):
        train_ds, val_ds, test_ds = self.input_split(smiles_list)
        self.save_data(train_ds, val_ds, test_ds)
        cmd = "cd {0}; bash fit.sh {1} {2} {3} {4} ; cd {5}".format(
            self.exec_folder,
            self.INPUT_DIR,
            self.OUTPUT_DIR,
            self.TRACE_DIR,
            self.CHECKPOINT_PATH,
            self.cwd,
        )
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)
        print("Fitting is done!")

    def _sample(self, n):
        print("Starting to Sample!")
        cmd = "cd {0}; python ./molecule_generation/cli/sample.py {1} {2} {3} ; cd {4}".format(
            self.exec_folder, self.tuned_checkpoint_dir, n, self.data_file, self.cwd
        )
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)
        print("Sampling is done!")

    def _read_molecules(self):
        output_file = self.data_file
        smiles_df = pd.read_csv(output_file)
        smiles_list = smiles_df["smiles"]
        return smiles_list

    def _clean(self):
        try:
            shutil.rmtree(self.tuned_checkpoint_dir)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

    def fit_and_sample(self, smiles_list, n):
        self._fit(smiles_list)
        self._sample(n)
        molecules = self._read_molecules()
        _clean()
        return molecules
