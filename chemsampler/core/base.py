import tempfile
import csv
import os
import pandas as pd

from ersilia import logger
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher
from ersilia.hub.delete.delete import ModelFullDeleter
from ersilia.hub.pull.pull import ModelPuller


class ModelArtifact(object):
    def __init__(self, model_id):
        self.model_id = model_id
        self.logger = logger
        try:
            self.load_model()
        except:
            self.model = None

    def load_model(self):
        self.model = ErsiliaModel(
            model=self.model_id,
            save_to_lake=False,
            service_class="pulled_docker",
            fetch_if_not_available=False,
        )

    def exists_locally(self):
        if self.model is None:
            return False
        else:
            return True

    def exists_remotely(self):
        mp = ModelPuller(model_id=self.model_id)
        return mp.is_available_in_dockerhub()

    def fetch(self):
        if not self.exists_locally():
            mf = ModelFetcher(force_from_dockerhub=True)
            mf.fetch(self.model_id)
            self.load_model()
        else:
            self.logger.info("Model {0} already exists".format(self.model_id))

    def delete(self):
        md = ModelFullDeleter()
        md.delete(self.model_id)

    def get_example_smiles_list(self):
        smiles_list = [
            "Nc1nc(NC2CC2)c2ncn([C@H]3C=C[C@@H](CO)C3)c2n1",
            "Nc1nc(=O)c2ncn(COCCO)c2[nH]1",
            "O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1",
            "NNC(=O)c1ccncc1",
            "CC/C(=C(\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
        ]
        return smiles_list

    def run(self, smiles_list):
        print(self.model_id)
        self.model.serve()
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        input_csv = os.path.join(tmp_folder, "input.csv")
        output_csv = os.path.join(tmp_folder, "output.csv")
        with open(input_csv, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles"])
            for smiles in smiles_list:
                writer.writerow([smiles])
        self.model.run(input=input_csv, output=output_csv)
        df = pd.read_csv(output_csv)
        print(df.head())
        self.model.close()
        return df
    
    def info(self):
        info = self.model.info()
        return info
