import os

from chemsampler.models.annotator import QEDAnnotator
from chemsampler.models.generator import HubGenerator
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import hill_climb

# Reserpine (PubChem CID 5770, InChIKey QEVHRUUCFGRFIF-MDEJGZGSSA-N):
# a complex, natural-product-derived antihypertensive drug. QED 0.374, MW 608.7.
SEED_SMILES = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

if __name__ == "__main__":
    summary, candidates_by_round = hill_climb(
        generator=HubGenerator("eos9taz"),
        annotators=[
            AnnotatorSpec("qed", QEDAnnotator(), cutoff=0.0, direction="higher")
        ],
        mode="sequential",
        seed_smiles=SEED_SMILES,
        n_rounds=5,
    )
    print(summary)

    os.makedirs("data", exist_ok=True)
    summary.to_csv("data/reserpine_qed_optimization_summary.csv", index=False)
    for round_num, candidates in candidates_by_round.items():
        candidates.to_csv(
            f"data/reserpine_qed_optimization_round{round_num}.csv", index=False
        )
