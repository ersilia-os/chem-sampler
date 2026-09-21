"""Build and locate the ChEMBL reference set used by :class:`~chemsampler.models.chembl.ChemblSampler`."""

import csv
import gzip
import hashlib
import os
import tempfile
import urllib.request
from multiprocessing import Pool

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from ..utils.logging import logger

RDLogger.DisableLog("rdApp.*")

CHEMBL_VERSION = "37"
CHEMREPS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    f"chembl_{CHEMBL_VERSION}_chemreps.txt.gz"
)
CHEMREPS_SHA256 = "ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7"

MW_MIN = 200.0
MW_MAX = 450.0
FILENAME = f"chembl_{CHEMBL_VERSION}_{MW_MIN:.0f}_{MW_MAX:.0f}Da.csv.gz"

_CHUNK = 20_000


def default_path() -> str:
    """
    Return the conventional location of the reference set.

    Returns
    -------
    str
        Path relative to the working directory, matching the layout `eosvc` expects.
    """
    return os.path.join("data", FILENAME)


def _filter_chunk(rows):
    """Keep single-component, RDKit-parseable rows whose molecular weight is in band."""
    kept = []
    for chembl_id, smiles in rows:
        if "." in smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        if MW_MIN <= mw <= MW_MAX:
            kept.append((chembl_id, smiles, f"{mw:.2f}"))
    return kept


def _chunks(path):
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        if header[:2] != ["chembl_id", "canonical_smiles"]:
            raise ValueError(f"Unexpected chemreps header: {header[:2]}")
        batch = []
        for row in reader:
            if len(row) < 2 or not row[1]:
                continue
            batch.append((row[0], row[1]))
            if len(batch) >= _CHUNK:
                yield batch
                batch = []
        if batch:
            yield batch


def _download(dest: str) -> None:
    logger.info(f"Downloading {CHEMREPS_URL}")
    urllib.request.urlretrieve(CHEMREPS_URL, dest)
    digest = hashlib.sha256()
    with open(dest, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    if digest.hexdigest() != CHEMREPS_SHA256:
        raise ValueError(
            f"Checksum mismatch for {dest}: expected {CHEMREPS_SHA256}, got {digest.hexdigest()}"
        )
    logger.success("Checksum verified")


def build(dest: str | None = None, n_processes: int | None = None) -> str:
    """
    Download ChEMBL and write the filtered reference set.

    Keeps single-component compounds whose RDKit molecular weight lies within
    [`MW_MIN`, `MW_MAX`], and writes them sorted by SMILES as a gzipped CSV with
    columns `chembl_id`, `smiles`, `mw`.

    Parameters
    ----------
    dest : str, optional
        Where to write the reference set. Defaults to :func:`default_path`.
    n_processes : int, optional
        Worker processes used for filtering. Defaults to one fewer than the CPU count.

    Returns
    -------
    str
        The path written.
    """
    dest = dest or default_path()
    n_processes = n_processes or max(1, (os.cpu_count() or 2) - 1)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="chemsampler-chembl-") as tmp_dir:
        raw = os.path.join(tmp_dir, "chemreps.txt.gz")
        _download(raw)

        logger.info(
            f"Filtering to {MW_MIN:.0f}-{MW_MAX:.0f} Da on {n_processes} processes"
        )
        kept = []
        with Pool(n_processes) as pool:
            for rows in pool.imap_unordered(_filter_chunk, _chunks(raw)):
                kept.extend(rows)

    kept.sort(key=lambda row: row[1])
    with gzip.open(dest, "wt", newline="", compresslevel=9) as f:
        writer = csv.writer(f)
        writer.writerow(["chembl_id", "smiles", "mw"])
        writer.writerows(kept)

    logger.success(f"Wrote {len(kept):,} compounds to {dest}")
    return dest


if __name__ == "__main__":
    build()
