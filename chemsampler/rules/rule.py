from rdkit import Chem
from rdkit.Chem import rdchem

class Ruler(object):

    def __init__(self,keep_smiles = None, avoid_smiles = None):
        self.keep_smiles = keep_smiles
        self.avoid_smiles = avoid_smiles

    def _convert_double_to_single_bonds(self,smiles):
        mol = Chem.MolFromSmiles(smiles)
        editable_mol = Chem.EditableMol(mol)
        for bond in mol.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                editable_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                editable_mol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), rdchem.BondType.SINGLE)
        return editable_mol.GetMol()

    def keep_substructure(self, sampled_smiles):
        keep_mol = self._convert_double_to_single_bonds(self.keep_smiles)
        sampled_smiles_filtered = []
        for smi in sampled_smiles:
            mol = self._convert_double_to_single_bonds(smi)
            if mol is not None:
                if mol.HasSubstructMatch(keep_mol):
                    sampled_smiles_filtered += [smi]
        return sampled_smiles_filtered

    def avoid_substructure(self, sampled_smiles):
        avoid_mol = self._convert_double_to_single_bonds(self.avoid_smiles)
        sampled_smiles_filtered = []
        for smi in sampled_smiles:
            mol = self._convert_double_to_single_bonds(smi)
            if mol is not None:
                if not mol.HasSubstructMatch(avoid_mol):
                    sampled_smiles_filtered += [smi]
        return sampled_smiles_filtered
