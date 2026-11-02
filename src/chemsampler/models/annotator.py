import pandas as pd

from .. import execution
from ..utils.logging import logger


class _RDKitAnnotator:
    """An in-house annotator: one RDKit score per molecule, from `_rdkit_scores`."""

    #: Its name in the annotator inputs (`source: rdkit`, `model: <name>`).
    name: str

    def score(self, smiles_list: list[str]) -> dict[str, float]:
        """
        Score a list of SMILES.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to score.

        Returns
        -------
        dict[str, float]
            Mapping from input SMILES to its score; SMILES RDKit cannot parse are
            left out.
        """
        ran = execution.current().run_rdkit(self.name, smiles_list)
        scores = {
            smi: float(value)
            for smi, value in zip(smiles_list, ran.table["value"])
            if pd.notna(value) and value != ""
        }
        logger.info(
            f"{self.name}: scored {len(smiles_list)} molecules ({ran.elapsed:.1f}s)"
        )
        return scores


class MolecularWeightAnnotator(_RDKitAnnotator):
    """
    Molecular weight scorer using RDKit: no Ersilia model needed.

    Scores are in Daltons (g/mol). Use with direction="higher" (cutoff sets minimum)
    or direction="lower" (cutoff sets maximum).

    Example: For drug-like MW range of 250-450 Da, use two annotators or set
    direction="higher" with cutoff=250 to enforce minimum.
    """

    name = "mw"


class QEDAnnotator(_RDKitAnnotator):
    """
    Drug-likeness scorer using RDKit's own QED implementation (Bickerton et al.,
    2012): no Ersilia model needed.

    Scores are on a 0-1 scale, higher is better.
    """

    name = "qed"


#: In-house annotators computed with RDKit, by the `model` name used in the
#: annotator inputs (`source: rdkit`).
RDKIT_ANNOTATORS = {"mw": MolecularWeightAnnotator, "qed": QEDAnnotator}
