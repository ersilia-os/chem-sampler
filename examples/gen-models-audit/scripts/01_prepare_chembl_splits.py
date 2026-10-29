"""
Build the seed set for the generative-models audit.

Downloads ChEMBL, keeps single-component compounds with 250 <= MW <= 450 Da,
draws a random sample of 1000, and splits it into 100 ten-compound CSV files
(one "smiles" column each) under results/ChEMBL_splits/. The raw download is
removed once the splits are written.

Fully self-contained: no imports from chemsampler, only stdlib + rdkit.
"""

import csv
import gzip
import os
import random
import shutil
import urllib.request

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

CHEMBL_VERSION = "37"
CHEMREPS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    f"chembl_{CHEMBL_VERSION}_chemreps.txt.gz"
)

MW_MIN = 250.0
MW_MAX = 450.0
N_COMPOUNDS = 1000
N_PER_FILE = 10
RANDOM_SEED = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
SPLITS_DIR = os.path.join(RESULTS_DIR, "ChEMBL_splits")
RAW_PATH = os.path.join(RESULTS_DIR, f"chembl_{CHEMBL_VERSION}_chemreps.txt.gz")


def download(dest: str) -> None:
    """Fetch the ChEMBL chemreps file to `dest`, creating parent dirs as needed."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print(f"Downloading {CHEMREPS_URL}")
    urllib.request.urlretrieve(CHEMREPS_URL, dest)


def filter_by_mw(path: str) -> list[str]:
    """Stream `path`, keeping single-component SMILES with MW in [MW_MIN, MW_MAX]."""
    kept = []
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        if header[:2] != ["chembl_id", "canonical_smiles"]:
            raise ValueError(f"Unexpected chemreps header: {header[:2]}")
        for row in reader:
            if len(row) < 2 or not row[1]:
                continue
            smiles = row[1]
            if "." in smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            if MW_MIN <= mw <= MW_MAX:
                kept.append(smiles)
    return kept


def write_splits(smiles: list[str]) -> None:
    """Write `smiles` as N_PER_FILE-row CSVs (column "smiles") into SPLITS_DIR."""
    if os.path.exists(SPLITS_DIR):
        shutil.rmtree(SPLITS_DIR)
    os.makedirs(SPLITS_DIR)

    n_files = len(smiles) // N_PER_FILE
    for i in range(n_files):
        chunk = smiles[i * N_PER_FILE : (i + 1) * N_PER_FILE]
        path = os.path.join(SPLITS_DIR, f"split_{i + 1:03d}.csv")
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles"])
            writer.writerows([[s] for s in chunk])
    print(f"Wrote {n_files} files of {N_PER_FILE} compounds to {SPLITS_DIR}")


def main():
    download(RAW_PATH)

    print(f"Filtering to single-component, {MW_MIN:.0f}-{MW_MAX:.0f} Da")
    filtered = filter_by_mw(RAW_PATH)
    print(f"{len(filtered):,} compounds in range")
    if len(filtered) < N_COMPOUNDS:
        raise ValueError(f"Only {len(filtered)} compounds in range, need {N_COMPOUNDS}")

    sample = random.Random(RANDOM_SEED).sample(filtered, N_COMPOUNDS)
    write_splits(sample)

    os.remove(RAW_PATH)
    print(f"Removed {RAW_PATH}")


if __name__ == "__main__":
    main()
