from preprocessing.main_preprocessor import preprocess_data
import pandas as pd
import os
import sys
from model.fine_tuner import FineTuner

root = os.path.dirname(os.path.abspath(__file__))


class _fit:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.data_file = os.path.join(root, "data")
        self.evaluation_file = os.path.join(root, "evaluation")

    def smiles_csv(self, smiles_list, filename):
        df = pd.DataFrame(smiles_list)
        df.to_csv(filename)

    def preprocess(self, smiles_list, name):
        filename_in_ = name + "_raw.csv"
        filename_in = os.path.join(self.data_file, filename_in_)
        filename_out_ = name + ".csv"
        filename_out = os.path.join(self.data_file, filename_out_)
        self.smiles_csv(smiles_list, filename_in)
        preprocess_data(
            filename_in=filename_in,
            filename_out=filename_out,
            model_type="ForwardRNN_512_fixed",
            starting_point="fixed",
            augmentation=1,
        )

    def fine_tune(self, model_name):
        for m in ["ForwardRNN_512_FineTuning_template"]:
            t = FineTuner(m)
            t.fine_tuning(model_name)

    # name_model is the file name for the fine tuned model
    def _fit(self, smiles_list, name_model):
        self.preprocess(smiles_list, name_model)
        self.fine_tune(name_model)


if __name__ == "__main__":
    smiles_list = sys.argv[1]
    name_model = sys.argv[2]
    _fit._fit(smiles_list, name_model)
