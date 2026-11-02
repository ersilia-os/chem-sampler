import pandas as pd

from .. import execution
from ..utils.logging import logger
from ..utils.similarity import tanimoto_to_seed


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


class TanimotoAnnotator:
    """
    Similarity to the run's seed: Morgan-fingerprint Tanimoto, on a 0-1 scale.

    Built without a seed; `hill_climb` binds the run's seed with `against`.
    Scored in the run's own process, never as a cluster job, so it needs no
    `rdkit` environment.

    Parameters
    ----------
    seed_smiles : str, optional
        The molecule candidates are compared with.
    """

    name = "tanimoto"

    def __init__(self, seed_smiles: str | None = None) -> None:
        self.seed_smiles = seed_smiles

    def against(self, seed_smiles: str) -> "TanimotoAnnotator":
        """A copy that compares candidates with `seed_smiles`."""
        return TanimotoAnnotator(seed_smiles)

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
            Mapping from input SMILES to its similarity to the seed; SMILES RDKit
            cannot parse are left out.

        Raises
        ------
        ValueError
            If no seed was bound.
        """
        if self.seed_smiles is None:
            raise ValueError(NEEDS_SEED)
        return tanimoto_to_seed(self.seed_smiles, smiles_list)


#: Why a Tanimoto annotator cannot be scored without a seed.
NEEDS_SEED = (
    "tanimoto compares candidates with the seed: give one "
    "(input/smiles.csv or --seed-smiles)"
)

#: In-house annotators computed with RDKit, by the `model` name used in the
#: annotator inputs (`source: rdkit`).
RDKIT_ANNOTATORS = {
    "mw": MolecularWeightAnnotator,
    "qed": QEDAnnotator,
    "tanimoto": TanimotoAnnotator,
}

#: Those of them a cluster job scores, with the `rdkit` environment.
RDKIT_ENV_ANNOTATORS = frozenset(
    name for name, cls in RDKIT_ANNOTATORS.items() if issubclass(cls, _RDKitAnnotator)
)
