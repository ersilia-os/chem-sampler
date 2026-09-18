import os

from chemsampler.models.annotator import QEDAnnotator
from chemsampler.models.generator import MolerGenerator
from chemsampler.optimize import hill_climb

# Reserpine (PubChem CID 5770, InChIKey QEVHRUUCFGRFIF-MDEJGZGSSA-N):
# a complex, natural-product-derived antihypertensive drug. QED 0.374, MW 608.7.
SEED_SMILES = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

if __name__ == "__main__":
    df = hill_climb(
        seed_smiles=SEED_SMILES,
        generator=MolerGenerator(),
        annotator=QEDAnnotator(),
        n_rounds=5,
    )
    print(df)
    os.makedirs("data", exist_ok=True)
    df.to_csv("data/reserpine_qed_optimization.csv", index=False)
