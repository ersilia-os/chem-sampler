from rdkit.Chem import Draw
from rdkit import Chem


mol_list = [
    "CCc1ccc(C(=O)N2CCCC2c2ccccc2)cc1-c1ccc(-c2ccccc2)cc1",
    "CCOC(=O)CCCCC(NC(=O)c1ccc(C)cc1)C(=O)N1CCCC1C(=O)NC(C(=O)C(F)(F)F)C(C)C",
    "CC#CCOc1ccc(S(=O)(=O)c2ccc3c(c2)N(CCN2CCNCC2)CC3)cc1",
]
m = Chem.MolFromSmiles(mol_list[0])
Draw.MolToImage(m)
