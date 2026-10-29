import pandas as pd
from rdkit import Chem


def candidates_from_output(df: pd.DataFrame, seed_smiles: str) -> list[str]:
    """
    Read the molecules out of a one-row Hub output table with `smi_*` columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of a Hub model run on the single compound `seed_smiles`.
    seed_smiles : str
        SMILES of the compound the model was run on.

    Returns
    -------
    list[str]
        Unique, RDKit-canonical SMILES in column order, without the seed itself
        and without anything RDKit cannot parse.
    """
    smi_columns = [c for c in df.columns if c.startswith("smi_")]
    raw_smiles = df.loc[0, smi_columns].dropna().tolist()

    seed_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(seed_smiles))
    candidates = {}
    for smi in raw_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical != seed_canonical:
            candidates[canonical] = None
    return list(candidates)
