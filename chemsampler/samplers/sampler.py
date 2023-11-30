import warnings

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import concurrent.futures

from ..core.base import ModelArtifact

warnings.filterwarnings("ignore")


class UnitSampler(ModelArtifact):
    def __init__(self, model_id, timeout_sec=1200):
        ModelArtifact.__init__(self, model_id=model_id)
        self.timeout_sec = timeout_sec

    def _calculate_similarity(self, ref_mol, mol_list):
        ref_fp = AllChem.GetMorganFingerprint(ref_mol, 2)
        mol_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in mol_list]
        similarities = [
            DataStructs.TanimotoSimilarity(ref_fp, mol_fp) for mol_fp in mol_fps
        ]
        return similarities

    def _sort_molecules_by_similarity(self, ref_mol, mol_list):
        similarities = self._calculate_similarity(ref_mol, mol_list)
        paired = list(zip(mol_list, similarities))
        sorted_mols = sorted(paired, key=lambda x: x[1], reverse=True)
        return [mol for mol, _ in sorted_mols]

    def _sort_by_similarity(self, smiles_list, ref_smiles):
        if len(smiles_list) < 2:
            return smiles_list
        ref_mol = Chem.MolFromSmiles(ref_smiles)
        mol_list = [Chem.MolFromSmiles(smi) for smi in smiles_list]
        sorted_mols = self._sort_molecules_by_similarity(
            ref_mol=ref_mol, mol_list=mol_list
        )
        return [Chem.MolToSmiles(mol) for mol in sorted_mols]

    def get_example_smiles(self):
        return self.get_example_smiles_list()[0]

    def sample(self, smiles):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(self.run, smiles_list=[smiles])
            try:
                df = future.result(timeout=self.timeout_sec)
                if df is not None:
                    for r in df.values:
                        r = r[2:]
                    sampled_smiles = []
                    for v in r:
                        try:
                            mol = Chem.MolFromSmiles(v)
                            smi = Chem.MolToSmiles(mol)
                        except:
                            smi = None
                        if smi is not None:
                            sampled_smiles += [smi]
                    sampled_smiles = list(set(sampled_smiles))
                    sampled_smiles = self._sort_by_similarity(sampled_smiles, smiles)
                    return sampled_smiles
                else:
                    self.logger.debug("Inner function timed out")
                    return []
            except concurrent.futures.TimeoutError:
                self.logger.debug("Outer function timed out")
                self.model.close()
                return []
            except:
                self.logger.debug("Some error occurred!")
                return []
    
    def get_info(self):
        info = self.info()
        return info