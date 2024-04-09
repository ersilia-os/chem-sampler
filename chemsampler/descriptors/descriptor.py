import numpy as np

from ..core.base import ModelArtifact


class DescriptorCalculator(ModelArtifact):
    def __init__(self, model_id):
        ModelArtifact.__init__(self, model_id=model_id)

    def calculate(self, smiles_list):
        df = self.run(smiles_list=smiles_list)
        X = np.array(df[list(df.columns)[2:]])
        return X
    
    def get_info(self):
        info = self.info()
        return info