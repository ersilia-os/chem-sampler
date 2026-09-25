import time

from rdkit import Chem

from ..hub.client import Backend, HubModel
from ..utils.logging import logger


class SeedRequiredError(ValueError):
    """Raised when a generator that needs a seed molecule is called without one."""


class HubGenerator:
    """
    Molecule generator backed by a generative model from the Ersilia Model Hub.

    Works with any Hub generative model that takes one compound and returns
    generated molecules in `smi_*` columns. How much of the input survives
    generation is model-specific: scaffold-based models rebuild from the input's
    Murcko scaffold, fragment-growing ones keep only a small ring fragment, and
    whole-molecule ones edit the input directly. Across the six default
    generators, mean Tanimoto to the input ranges from about 0.13 to 0.68
    (benchmark in ersilia-os/ersilia#1919).

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the generative model (e.g. "eos9taz").
    backend : {"ersilia", "run_sh"}, optional
        Passed through to `HubModel`, by default "run_sh".
    """

    def __init__(self, model_id: str, backend: Backend = "run_sh"):
        self.model_id = model_id
        self._hub_model = HubModel(model_id, backend=backend)

    def generate(self, seed_smiles: str | None) -> list[str]:
        """
        Generate candidate molecules from a seed molecule.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        list[str]
            Unique, RDKit-canonicalized candidate SMILES, excluding the seed itself
            and any outputs RDKit cannot parse.

        Raises
        ------
        SeedRequiredError
            If `seed_smiles` is `None`. Hub generative models take one compound as
            input, so there is no seedless mode for them.
        """
        if seed_smiles is None:
            raise SeedRequiredError(f"{self.model_id} requires a seed molecule")

        df = self._hub_model.run([seed_smiles])
        smi_columns = [c for c in df.columns if c.startswith("smi_")]
        raw_smiles = df.loc[0, smi_columns].dropna().tolist()

        seed_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(seed_smiles))
        candidates = set()
        for smi in raw_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canonical = Chem.MolToSmiles(mol)
            if canonical != seed_canonical:
                candidates.add(canonical)
        return list(candidates)

    def generate_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        """
        Generate candidates, keyed by this generator's model id.

        Gives `HubGenerator` the same provenance-tracking interface as
        `GeneratorPool`, so callers can treat a single generator and a pool alike.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        dict[str, list[str]]
            `{self.model_id: self.generate(seed_smiles)}`.
        """
        return {self.model_id: self.generate(seed_smiles)}


class GeneratorPool:
    """
    Several generators used as one, pooling their candidates.

    Satisfies the same interface as `HubGenerator`, so it can be passed wherever
    a single generator is expected.

    Parameters
    ----------
    generators : list[HubGenerator]
        Generators whose candidates are pooled.
    """

    def __init__(self, generators: list[HubGenerator]):
        self.generators = generators

    def generate(self, seed_smiles: str | None) -> list[str]:
        """
        Generate candidates from every generator and return their union.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        list[str]
            Deduplicated union of all generators' candidates.
        """
        candidates = set()
        for generator in self.generators:
            result = self._generate_one(generator, seed_smiles)
            candidates.update(result)
        return list(candidates)

    def generate_by_model(self, seed_smiles: str | None) -> dict[str, list[str]]:
        """
        Generate candidates, keeping track of which model produced each set.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        dict[str, list[str]]
            Mapping from model identifier to that model's candidates.
        """
        by_model = {}
        for generator in self.generators:
            by_model[generator.model_id] = self._generate_one(generator, seed_smiles)
        return by_model

    @staticmethod
    def _generate_one(generator, seed_smiles: str | None) -> list[str]:
        """Run one generator, turning a missing-seed error into a warning + skip."""
        model_id = getattr(generator, "model_id", generator)
        start = time.monotonic()
        try:
            result = generator.generate(seed_smiles)
        except SeedRequiredError:
            logger.warning(f"{model_id}: requires a seed molecule, skipped")
            return []
        elapsed = time.monotonic() - start
        if not result:
            logger.warning(
                f"{model_id}: contributed 0 candidates for this seed ({elapsed:.1f}s)"
            )
        else:
            logger.info(
                f"{model_id}: contributed {len(result)} candidates for this seed "
                f"({elapsed:.1f}s)"
            )
        return result
