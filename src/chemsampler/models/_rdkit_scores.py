"""The in-house RDKit annotators' scores, in a file that needs only RDKit.

chemsampler runs a copy of this file with the `rdkit` environment's Python
(the one `chemsampler fetch rdkit` builds), here or as a cluster job:

    python _rdkit_scores.py <name> <input.csv> <output.csv> [<seed>]

The input has a `smiles` column; the output has `smiles` and `value`, empty
where RDKit cannot parse the SMILES. `tanimoto` needs the seed. Keep this file
free of chemsampler imports.
"""

import csv
import sys

from rdkit import Chem
from rdkit.Chem import DataStructs, Descriptors, rdFingerprintGenerator

#: Radius-2, 2048-bit Morgan fingerprints (ECFP4-like), for `tanimoto`.
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

#: Each annotator's name in the annotator inputs (`source: rdkit`).
NAMES = ("mw", "qed", "tanimoto")


def score(name: str, smiles_list: list[str], seed: str | None = None) -> list:
    """One value per SMILES, in order; None where RDKit cannot parse it."""
    if name == "tanimoto":
        return tanimoto(seed, smiles_list)
    scorer = {"mw": Descriptors.MolWt, "qed": Descriptors.qed}[name]
    values = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        values.append(None if mol is None else scorer(mol))
    return values


def tanimoto(seed: str | None, smiles_list: list[str]) -> list:
    """Each SMILES' Morgan-fingerprint Tanimoto similarity to `seed`, 0 to 1."""
    seed_mol = None if seed is None else Chem.MolFromSmiles(seed)
    if seed_mol is None:
        raise ValueError(f"tanimoto needs a seed RDKit can read, got {seed!r}")
    seed_fp = _MORGAN.GetFingerprint(seed_mol)
    values = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        values.append(
            None
            if mol is None
            else DataStructs.TanimotoSimilarity(seed_fp, _MORGAN.GetFingerprint(mol))
        )
    return values


def read(path) -> list[float | None]:
    """The values `main` wrote, exactly (pandas' default float parser is not)."""
    with open(path, newline="") as f:
        return [
            float(row["value"]) if row["value"] else None for row in csv.DictReader(f)
        ]


def main(argv: list[str]) -> None:
    name, source, target, *seed = argv
    with open(source, newline="") as f:
        smiles_list = [row["smiles"] for row in csv.DictReader(f)]
    values = score(name, smiles_list, *seed)
    with open(target, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "value"])
        for smiles, value in zip(smiles_list, values):
            writer.writerow([smiles, "" if value is None else repr(value)])


if __name__ == "__main__":
    main(sys.argv[1:])
