from ...tools.bimodal.sampler import _BimodalSampler
from ...tools.fpsim2.searcher import SimilaritySearcher, RandomSearcher
import os

root = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.abspath(
    os.path.join(root, "..", "..", "..", "data", "pre_calculated", "bimodal")
)


class BimodalSampler:
    def __init__(self):
        # TODO: Come up with a better way to store precalculated files
        self.fp_filename = os.path.join(data_path, "bimodal_10k.h5")
        self.db_smiles_filename = os.path.join(data_path, "bimodal_10k.csv")

        if not os.path.exists(self.fp_filename):
            smiles_list = SimilaritySearcher(self.fp_filename).read_db_smiles()
            SimilaritySearcher(self.fp_filename).fit(smiles_list)

    def sample(self, smiles_list, n, search_pre_calculated=False, cutoff=0.4):
        if search_pre_calculated:
            samples = []
            for smile in smiles_list:
                _samples = SimilaritySearcher(self.fp_filename).search(smile, cutoff)
                if len(_samples) == 0:
                    _samples = []
                samples += _samples

        else:
            sampler = _BimodalSampler()
            samples = sampler.sample(n)

        return samples
