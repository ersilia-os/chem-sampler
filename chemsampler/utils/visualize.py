from rdkit import Chem
from rdkit.Chem import Draw
import os
from .config import ConfigRun

class VisualizeMolecules(object):
    def __init__(self, num_samples =100):
        self.num_samples= num_samples

    def visualize_mols(self, df, output_folder, round):
        max = df['sampled_smiles'].nunique()
        if max < self.num_samples:
            sampled_smiles = df['sampled_smiles'].sample(n=max).tolist()
        else:
            sampled_smiles = df['sampled_smiles'].sample(n=self.num_samples).tolist()
        sampled_mols = [Chem.MolFromSmiles(smile) for smile in sampled_smiles]
        img = Draw.MolsToGridImage(mols=sampled_mols, molsPerRow=5)
        img.save(os.path.join(output_folder, "mols_round{}.png".format(round)))
