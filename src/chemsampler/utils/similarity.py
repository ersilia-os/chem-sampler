from ..models import _rdkit_scores


def tanimoto_to_seed(seed_smiles: str, smiles_list: list[str]) -> dict[str, float]:
    """
    Morgan-fingerprint Tanimoto similarity of each candidate to a seed molecule.

    Uses radius-2, 2048-bit fingerprints (ECFP4-equivalent): the same values as
    the `tanimoto` annotator, computed here with chemsampler's own RDKit.

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
    values = _rdkit_scores.tanimoto(seed_smiles, smiles_list)
    return {smi: v for smi, v in zip(smiles_list, values) if v is not None}
