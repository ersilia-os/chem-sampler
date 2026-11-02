"""The in-house RDKit annotators' scores, in a file that needs only RDKit.

chemsampler scores with `score`, and a cluster job runs this same file with any
Python that has RDKit (the `rdkit` environment `chemsampler fetch` builds), so
both give the same values:

    python _rdkit_scores.py <name> <input.csv> <output.csv>

The input has a `smiles` column; the output has `smiles` and `value`, empty
where RDKit cannot parse the SMILES. Keep this file free of chemsampler imports.
"""

import csv
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

#: Each annotator's name in the annotator inputs (`source: rdkit`) and its score.
SCORERS = {"mw": Descriptors.MolWt, "qed": Descriptors.qed}


def score(name: str, smiles_list: list[str]) -> list[float | None]:
    """One value per SMILES, in order; None where RDKit cannot parse it."""
    scorer = SCORERS[name]
    values = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        values.append(None if mol is None else scorer(mol))
    return values


def read(path) -> list[float | None]:
    """The values `main` wrote, exactly (pandas' default float parser is not)."""
    with open(path, newline="") as f:
        return [
            float(row["value"]) if row["value"] else None for row in csv.DictReader(f)
        ]


def main(argv: list[str]) -> None:
    name, source, target = argv
    with open(source, newline="") as f:
        smiles_list = [row["smiles"] for row in csv.DictReader(f)]
    with open(target, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "value"])
        for smiles, value in zip(smiles_list, score(name, smiles_list)):
            writer.writerow([smiles, "" if value is None else repr(value)])


if __name__ == "__main__":
    main(sys.argv[1:])
