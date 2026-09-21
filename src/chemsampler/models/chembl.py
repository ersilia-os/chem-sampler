import csv
import gzip
import os
import random
import shutil
import subprocess

from ..data.chembl import default_path
from ..utils.logging import logger


class ChemblSampler:
    """
    Null baseline that draws molecules at random from a ChEMBL reference set.

    Satisfies the same interface as `HubGenerator`, so it can be passed wherever a
    generator is expected. It is **not** a generative model: `generate` ignores the
    seed entirely and returns an unconditioned random sample. It is therefore
    deliberately excluded from `VALIDATED_GENERATORS`, and must be opted into
    explicitly.

    It exists to answer "does a generator beat drawing at random from known
    chemistry?". On seed-agnostic objectives such as QED a random draw is hard to
    beat, precisely because it is unconstrained by the seed.

    Parameters
    ----------
    n : int, optional
        Number of molecules to draw per call, by default 1000.
    random_state : int, optional
        Seed for the sampling RNG. If None, draws differ between calls.
    path : str, optional
        Location of the reference set. Defaults to
        :func:`chemsampler.data.chembl.default_path`. If the file is missing it is
        fetched with `eosvc`.
    """

    model_id = "chembl"

    def __init__(
        self, n: int = 1000, random_state: int | None = None, path: str | None = None
    ):
        self.n = n
        self.random_state = random_state
        self.path = path or default_path()
        self._smiles = None

    def _fetch(self) -> None:
        """Pull the reference set from S3 with the eosvc CLI."""
        if shutil.which("eosvc") is None:
            raise FileNotFoundError(
                f"ChEMBL reference set not found at {self.path} and the eosvc CLI is "
                "not available. Install eosvc, or build the set locally with "
                "`python -m chemsampler.data.chembl`."
            )
        logger.info(f"Reference set missing, fetching with eosvc: {self.path}")
        result = subprocess.run(
            ["eosvc", "download", "--path", self.path],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0 or not os.path.exists(self.path):
            raise RuntimeError(
                f"eosvc could not fetch {self.path}: {result.stderr.strip() or 'unknown error'}"
            )

    def _load(self) -> list[str]:
        """Read and cache the SMILES column, fetching the file if it is absent."""
        if self._smiles is None:
            if not os.path.exists(self.path):
                self._fetch()
            with gzip.open(self.path, "rt", newline="") as f:
                reader = csv.reader(f)
                next(reader)
                self._smiles = [row[1] for row in reader]
            logger.info(f"Loaded {len(self._smiles):,} ChEMBL compounds")
        return self._smiles

    def generate(self, seed_smiles: str) -> list[str]:
        """
        Draw `n` molecules at random, ignoring the seed.

        Parameters
        ----------
        seed_smiles : str
            Accepted for interface compatibility and not used.

        Returns
        -------
        list[str]
            `n` SMILES sampled without replacement from the reference set.
        """
        smiles = self._load()
        if self.n > len(smiles):
            raise ValueError(
                f"Requested {self.n} molecules but the reference set holds {len(smiles)}."
            )
        return random.Random(self.random_state).sample(smiles, self.n)
