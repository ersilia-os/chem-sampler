from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator

_MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def tanimoto_to_seed(seed_smiles: str, smiles_list: list[str]) -> dict[str, float]:
    """
    Morgan-fingerprint Tanimoto similarity of each candidate to a seed molecule.

    Uses radius-2, 2048-bit fingerprints (ECFP4-equivalent).

    Parameters
    ----------
    seed_smiles : str
        SMILES of the seed molecule.
    smiles_list : list[str]
        Candidate SMILES to compare against the seed.

    Returns
    -------
    dict[str, float]
        Mapping from candidate SMILES to its Tanimoto similarity to the seed, on a
        0-1 scale. Candidates RDKit cannot parse are omitted.

    Raises
    ------
    ValueError
        If `seed_smiles` cannot be parsed by RDKit.
    """
    seed_mol = Chem.MolFromSmiles(seed_smiles)
    if seed_mol is None:
        raise ValueError(f"Could not parse seed SMILES: {seed_smiles!r}")
    seed_fp = _MORGAN_GENERATOR.GetFingerprint(seed_mol)

    similarities = {}
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = _MORGAN_GENERATOR.GetFingerprint(mol)
        similarities[smi] = DataStructs.TanimotoSimilarity(seed_fp, fp)
    return similarities
