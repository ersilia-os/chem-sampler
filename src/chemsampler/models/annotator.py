from rdkit import Chem
from rdkit.Chem import Descriptors


class QEDAnnotator:
    """
    Drug-likeness scorer using RDKit's own QED implementation (Bickerton et al.,
    2012). Computed locally - no Ersilia model or Docker/network dependency.

    Scores are on a 0-1 scale, higher is better.
    """

    def score(self, smiles_list: list[str]) -> dict[str, float]:
        """
        Score a list of SMILES by QED.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to score.

        Returns
        -------
        dict[str, float]
            Mapping from input SMILES to its QED score.
        """
        scores = {}
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                scores[smi] = Descriptors.qed(mol)
        return scores
