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
        self.CHECKPOINT_PATH = os.path.join(root, "MODEL_DIR","GNN_Edge_MLP_MoLeR__2022-02-24_07-16-23_best.pkl"  )
        self.INPUT_DIR = tempfile.mkdtemp()
        self.OUTPUT_DIR = tempfile.mkdtemp()
        self. TRACE_DIR = tempfile.mkdtemp()
        self.train_split = 0.8
        self.test_split = 0.1
        self.val_split = 0.1 
    
    def input_split(self, smiles_list):
        assert (self.train_split + self.test_split + self.val_split) == 1
        assert self.val_split == self.test_split 
        df = pd.DataFrame(smiles_list)
        df_sample = df.sample(frac=1, random_state=12)
        indices_or_sections = [int(self.train_split * len(df)), int((1 - self.val_split - self.test_split) * len(df))]
        train_ds, val_ds, test_ds = np.split(df_sample, indices_or_sections)
        return train_ds, val_ds, test_ds

    def write_smi(self, smiles, filename):
        csv_filename = filename + '.csv'
        csv_filename = os.path.join(self.INPUT_DIR, csv_filename)
        smi_filename = filename + '.smiles'
        smi_filename = os.path.join(self.INPUT_DIR, smi_filename)
        smiles.to_csv(csv_filename, header=False, index=False)
        os.rename(csv_filename, smi_filename)
      
    def save_data(self, train_ds, val_ds, test_ds):
        self.write_smi(train_ds, 'train')
        self.write_smi(val_ds, 'valid')
        self.write_smi(test_ds, 'test') 
        
    def fit(self, smiles_list):
        train_ds, val_ds, test_ds = self.input_split(smiles_list)
        self.save_data(train_ds, val_ds, test_ds)
        cmd = "cd {0}; bash fit.sh {1} {2} {3} {4} ; cd {5}".format(
            self.exec_folder, self.INPUT_DIR , self.OUTPUT_DIR , self.TRACE_DIR, self.CHECKPOINT_PATH, self.cwd
        )
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)




