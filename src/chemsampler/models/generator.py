from rdkit import Chem

from ..hub.client import HubModel


class MolerGenerator:
    """
    Scaffold-preserving molecule generator backed by eos9taz (MoLeR + Enamine fragments).

    Given a seed molecule, generates up to 1000 unique candidates that extend its
    scaffold with fragments sampled from Enamine's library. The candidate count is
    not configurable - it is intrinsic to the model's output width.
    """

    MODEL_ID = "eos9taz"

    def __init__(self):
        self._hub_model = HubModel(self.MODEL_ID)

    def generate(self, seed_smiles: str) -> list[str]:
        """
        Generate scaffold-preserving candidates for a seed molecule.

        Parameters
        ----------
        seed_smiles : str
            SMILES of the seed molecule.

        Returns
        -------
        list[str]
            Unique, RDKit-canonicalized candidate SMILES, excluding the seed itself
            and any invalid outputs.
        """
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
