from ..hub.client import SIMILARITY_MODEL_IDS, Backend, HubModel
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
        Ersilia identifier of the similarity model. Only "eos8nyg" (Ersilia
        Reference Library Similarity) is supported for now.
    backend : {"ersilia", "run_sh"}, optional
        Passed through to `HubModel`, by default "run_sh".

    Raises
    ------
    ValueError
        If `model_id` is not a supported similarity model.
    """

    def __init__(self, model_id: str, backend: Backend = "run_sh"):
        if model_id not in SIMILARITY_MODEL_IDS:
            supported = ", ".join(sorted(SIMILARITY_MODEL_IDS))
            raise ValueError(
                f"{model_id} is not a supported similarity model (supported: {supported})"
            )
        self.model_id = model_id
        self._hub_model = HubModel(model_id, backend=backend)

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

        df = self._hub_model.run([seed_smiles])
        return candidates_from_output(df, seed_smiles)

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
