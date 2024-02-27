import pandas as pd

class InputSelector(object):
    def __init__(self, info_file, results_file, saturation_number):
        self.info = info_file
        self.results = results_file
        self.saturation_number = saturation_number
    
    def _is_not_saturated(self):
        new_cpds = int(self.info[-1]["total_unique_generated"])
        print(new_cpds)
        if new_cpds > self.saturation_number:
            return True
        else:
            return False
        
    def _rank_columns(self, df):
        euclidean_cols = df.filter(like='_euclidean')
        tanimoto_cols = df.filter(like='_tanimoto')
        euclidean_rank_cols = euclidean_cols.rank(axis=1)
        tanimoto_rank_cols = tanimoto_cols.rank(axis=1, ascending=False)
        all_rank_columns = pd.concat([euclidean_rank_cols, tanimoto_rank_cols], axis=1)
        df['avg_rank'] = all_rank_columns.mean(axis=1)
        best_smiles = df.loc[df['avg_rank'].idxmin(), 'sampled_smiles']
        print(best_smiles)
        return best_smiles


    def choose_input(self):
        if len(self.info)==1:
            input_smiles = self.info[0]["seed_smiles"]
        else:
            if self._is_not_saturated():
                input_smiles = self.info[-1]["input_smiles"]
            else:
                results = self.results[self.results["round"]==self.info[-1]["round_number"]]
                input_smiles = self._rank_columns(results)
        return input_smiles