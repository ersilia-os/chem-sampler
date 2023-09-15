from ...tools.mgm.sampler import _MgmSampler
from ...tools.fpsim2.searcher import SimilaritySearcher, RandomSearcher
import os


class MgmSampler:
    def __init__(self):
        # #TODO: Come up with a better way to store precalculated files
        # self.fp_filename = '../data/pre_calculated/moler/moler_100k.h5'
        # self.db_smiles_filename = '../data/pre_calculated/moler/moler_10k.csv'

        # if os.path.exists(self.fp_filename) is False:
        #     smiles_list = SimilaritySearcher(self.fp_filename).read_db_smiles()
        #     SimilaritySearcher(self.fp_filename).fit(smiles_list)
        pass

    def sample(self, n, smiles_list=[], search_pre_calculated=False, cutoff=0.7):
        if search_pre_calculated == True:
            samples = []
            for smile in smiles_list:
                _samples = SimilaritySearcher(self.fp_filename).search(smile, cutoff)
                if len(_samples) == 0:
                    print(
                        "Coudn't find samples in the cutoff range : {}. Generating random samples:".format(
                            cutoff
                        )
                    )
                    _samples = RandomSearcher(self.db_smiles_filename).search(n)
                samples += _samples

        else:
            sampler = _MgmSampler()
            samples = sampler.sample(n)

        return samples
