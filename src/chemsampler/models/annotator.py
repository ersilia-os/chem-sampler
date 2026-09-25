from rdkit import Chem
from rdkit.Chem import Descriptors


class MolecularWeightAnnotator:
    """
    Molecular weight scorer using RDKit. Computed locally - no Ersilia model or
    network dependency.

    Scores are in Daltons (g/mol). Use with direction="higher" (cutoff sets minimum)
    or direction="lower" (cutoff sets maximum).

    Example: For drug-like MW range of 250-450 Da, use two annotators or set
    direction="higher" with cutoff=250 to enforce minimum.
    """

    def score(self, smiles_list: list[str]) -> dict[str, float]:
        """
        Score a list of SMILES by molecular weight.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to score.

        Returns
        -------
        dict[str, float]
            Mapping from input SMILES to its molecular weight in Da.
        """
        scores = {}
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                scores[smi] = Descriptors.MolWt(mol)
        return scores


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
