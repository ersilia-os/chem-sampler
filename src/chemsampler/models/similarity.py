import time

from ..hub.client import (
    SIMILARITY_MODEL_IDS,
    HubModel,
    not_allowed_similarity_message,
)
from ..utils.logging import logger
from ._smiles import candidates_from_output
from .generator import SeedRequiredError


class HubSimilarity:
    """
    Candidate source backed by a similarity-search model from the Ersilia Model Hub.

    Looks the seed up in a reference library and returns its nearest neighbours.
    These are known compounds, not newly generated ones, and the search is
    deterministic. It is independent of `HubGenerator`: it is not a generator and
    has no `generate_by_model`, so it cannot be passed where a generator is
    expected. Pass it to `hill_climb(similarity=...)` instead.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the similarity model: "eos8nyg" (Ersilia Reference
        Library), "eos1bgt" (COCONUT natural products) or "eos80al" (Enamine Hit
        Locator), see `SIMILARITY_MODEL_IDS`.
    models_dir : str
        Folder holding the model repositories, passed through to `HubModel`.
    envs_dir : str, optional
        Folder holding the models' conda environments, passed through to
        `HubModel`.

    Raises
    ------
    ValueError
        If `model_id` is not a supported similarity model.
    """

    def __init__(
        self,
        model_id: str,
        models_dir: str | None = None,
        envs_dir: str | None = None,
    ):
        if model_id not in SIMILARITY_MODEL_IDS:
            raise ValueError(not_allowed_similarity_message(model_id))
        self.model_id = model_id
        self._hub_model = HubModel(model_id, models_dir=models_dir, envs_dir=envs_dir)

    def search(self, seed_smiles: str | None) -> list[str]:
        """
        Find the seed's nearest neighbours in the model's reference library.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        list[str]
            Unique, RDKit-canonical SMILES, nearest first, excluding the seed
            itself and any outputs RDKit cannot parse.

        Raises
        ------
        SeedRequiredError
            If `seed_smiles` is `None`: a similarity search needs a query molecule.
        """
        if seed_smiles is None:
            raise SeedRequiredError(f"{self.model_id} requires a seed molecule")

        start = time.monotonic()
        df = self._hub_model.run([seed_smiles])
        candidates = candidates_from_output(df, seed_smiles)
        found = (
            f"{self.model_id}: found {len(candidates)} neighbours for this seed "
            f"({time.monotonic() - start:.1f}s)"
        )
        if candidates:
            logger.info(found)
        else:
            logger.warning(found)
        return candidates

    def search_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        """
        Search, keyed by this model's id so candidates keep their provenance.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        dict[str, list[str]]
            `{self.model_id: self.search(seed_smiles)}`.
        """
        return {self.model_id: self.search(seed_smiles)}
