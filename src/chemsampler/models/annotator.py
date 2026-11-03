import os

import pandas as pd

from .. import execution
from ..hub.envs import RDKIT_ENV, env_path, is_complete
from ..utils.logging import logger


class _RDKitAnnotator:
    """
    An in-house annotator: one RDKit score per molecule, from `_rdkit_scores`.

    It runs in the `rdkit` environment (`chemsampler fetch rdkit`), here as its
    own process or on a cluster as a job, like a Hub model in its own.

    Parameters
    ----------
    envs_dir : str, optional
        The folder of environments that holds `rdkit`. Without it (an annotator
        built by hand in Python), the score is computed with chemsampler's own
        RDKit, which gives the same values.
    """

    #: Its name in the annotator inputs (`source: rdkit`, `model: <name>`).
    name: str

    def __init__(self, envs_dir: str | None = None) -> None:
        self.envs_dir = envs_dir

    def env_dir(self):
        """
        The `rdkit` environment, or None without an environments folder.

        Raises
        ------
        RuntimeError
            If an environments folder is set but its `rdkit` is not complete.
        """
        if not self.envs_dir:
            return None
        env_dir = env_path(os.path.expanduser(self.envs_dir), RDKIT_ENV)
        if not is_complete(env_dir):
            raise RuntimeError(
                f"{self.name}: no complete rdkit environment at {env_dir}; run "
                "`chemsampler fetch rdkit`"
            )
        return env_dir

    def args(self) -> list[str]:
        """What `_rdkit_scores` needs besides the molecules (the seed, for one)."""
        return []

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
        ran = execution.current().run_rdkit(self, smiles_list)
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


class TanimotoAnnotator(_RDKitAnnotator):
    """
    Similarity to the run's seed: Morgan-fingerprint Tanimoto, on a 0-1 scale.

    Built without a seed; `hill_climb` binds the run's seed with `against`.

    Parameters
    ----------
    seed_smiles : str, optional
        The molecule candidates are compared with.
    envs_dir : str, optional
        As for every RDKit annotator.
    """

    name = "tanimoto"

    def __init__(
        self, seed_smiles: str | None = None, envs_dir: str | None = None
    ) -> None:
        super().__init__(envs_dir)
        self.seed_smiles = seed_smiles

    def against(self, seed_smiles: str) -> "TanimotoAnnotator":
        """A copy that compares candidates with `seed_smiles`."""
        return TanimotoAnnotator(seed_smiles, self.envs_dir)

    def args(self) -> list[str]:
        if self.seed_smiles is None:
            raise ValueError(NEEDS_SEED)
        return [self.seed_smiles]


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
