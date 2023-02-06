from ...tools.fast_jtnn.sampler import _JtnnSampler

from ...tools.fpsim2.searcher import SimilaritySearcher, RandomSearcher
import os

class JtnnSampler:
    def __init__(self):
        #TODO: pre-calculate samples
        self.fp_filename = '../chemsampler/data/pre_calculated/fast_jtnn/.h5'
        self.db_smiles_filename = '../chemsampler/data/pre_calculated/fast_jtnn/.csv'

        if os.path.exists(self.fp_filename) is False:
            smiles_list = SimilaritySearcher(self.fp_filename).read_db_smiles()
            SimilaritySearcher(self.fp_filename).fit(smiles_list)

    def sample(self, smiles_list, n,  search_pre_calculated=True, cutoff = 0.7):
        if search_pre_calculated==True:
            samples = []
            for smile in smiles_list:
                _samples = SimilaritySearcher(self.fp_filename).search(smile, cutoff=0.7)
                if len(_samples) == 0:
                    print('Coudn\'t find samples in the cutoff range : {}. Generating random samples:'.format(cutoff))
                    _samples = RandomSearcher(self.db_smiles_filename).search(n)
                samples += _samples 
        
        else:
            sampler = _JtnnSampler()
            samples = sampler.sample(n)

        return samples


