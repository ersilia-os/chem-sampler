import os
import csv
from .sampler import ChemSampler


root = os.path.abspath(os.path.dirname(__file__))


def example():
    smiles_list = []
    with open(os.path.join(root,"..", "data", "example_drugs.csv"), "r") as f:
        reader = csv.reader(f)
        next(reader)
        for r in reader:
            smiles_list += [r[0]]
    return smiles_list
