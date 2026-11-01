import time

from rdkit import Chem
from rdkit.Chem import Descriptors

from ..utils.logging import logger


def _report(name: str, count: int, start: float) -> None:
    """Log how many molecules an annotator scored and how long it took."""
    logger.info(f"{name}: scored {count} molecules ({time.monotonic() - start:.1f}s)")


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
        start = time.monotonic()
        scores = {}
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                scores[smi] = Descriptors.MolWt(mol)
        _report("mw", len(smiles_list), start)
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
        start = time.monotonic()
        scores = {}
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                scores[smi] = Descriptors.qed(mol)
        _report("qed", len(smiles_list), start)
        return scores


#: In-house annotators computed with RDKit, by the `model` name used in the
#: annotator inputs (`source: rdkit`).
RDKIT_ANNOTATORS = {"mw": MolecularWeightAnnotator, "qed": QEDAnnotator}
