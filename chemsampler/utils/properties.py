import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

class PropertyCalculator(object):
    def __init__(self, seed_smiles, sampled_smiles):
        self.seed_smiles = seed_smiles
        self.sampled_smiles = sampled_smiles

    def calc_molecular_weight(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.sampled_smiles]
        weights = [Descriptors.MolWt(mol) if mol is not None else None for mol in mols]
        return weights

    def calc_qed(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.sampled_smiles]
        qed_values = [Descriptors.qed(mol) if mol is not None else None for mol in mols]
        return qed_values

    def calc_logp(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.sampled_smiles]
        logp_values = [Descriptors.MolLogP(mol) if mol is not None else None for mol in mols]
        return logp_values

    def run(self):
        mw = self.calc_molecular_weight()
        qed = self.calc_qed()
        logp = self.calc_logp()
        df = pd.DataFrame({
            "sampled_smiles": self.sampled_smiles,
            "mw": mw,
            "qed": qed,
            "logp": logp
        })
        return df
