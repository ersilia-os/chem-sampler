from ...utils.conda import SimpleConda
import os

root = os.path.dirname(os.path.abspath(__file__))

BIMODAL_CONDA_ENVIRONMENT = "bimodal"


class fit:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)

    def fit(self, smiles_list, out_file_name):
        cmd = "cd {0} ; python _fit.py {1} {2}; cd {3}".format(
            self.exec_folder, smiles_list, out_file_name, self.cwd
        )
        SimpleConda().run_commandlines(BIMODAL_CONDA_ENVIRONMENT, cmd)
