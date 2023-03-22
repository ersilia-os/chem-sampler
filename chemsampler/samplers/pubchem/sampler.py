import requests
import random
from tqdm import tqdm
from timeit import default_timer as timer
from ratelimit import sleep_and_retry, limits
from standardiser import standardise
from rdkit import Chem


@sleep_and_retry
def run_chemed(origin_smiles: str, num_samples: int, similarity: float = 0.7):
    """Function adapted from Andrew White's Exmol"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{requests.utils.quote(origin_smiles)}/property/CanonicalSMILES/JSON"
    try:
        reply = requests.get(
            url,
            params={"Threshold": int(similarity *100), "MaxRecords": num_samples},
            headers={"accept": "text/json"},
            timeout=10,
        )
    except requests.exceptions.Timeout:
        print("Pubchem seems to be down right now")
        return []
    try:
        data = reply.json()
        if "PropertyTable" not in data or "Properties" not in data["PropertyTable"]:
            return []
        smiles = [d["CanonicalSMILES"] for d in data["PropertyTable"]["Properties"]]
        smiles = list(set(smiles))
        return smiles
    except:
        return []
    


class PubChemSampler(object):
    def __init__(self):
        self.inflation = 10

    def _sample(self, smiles, n):
        smiles = run_chemed(origin_smiles=smiles, num_samples=n)
        std_smiles = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def sample(self, smiles_list, n, time_budget_sec=600):
        n_individual = int(n / len(smiles_list))
        n_individual = int(n_individual * self.inflation)
        sampled_smiles = []
        t0 = timer()
        for smi in tqdm(smiles_list):
            sampled_smiles += self._sample(smi, n_individual)
            t1 = timer()
            if (t1 - t0) > time_budget_sec:
                break
        sampled_smiles = list(set(sampled_smiles))
        random.shuffle(sampled_smiles)
        return sampled_smiles


