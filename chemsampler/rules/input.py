from ..utils.config import ConfigRun

class InputSelector(object):
    def __init__(self, config_file):
        rc = ConfigRun(config_file)
        self.info = rc.load_info_file()
        self.results = rc.load_results()
    
    def _is_not_saturated(self):
        new_cpds = int(self.info[-1]["unique_generated"])
        if new_cpds > 200:
            return True
        else:
            return False
        
    def choose_input(self):
        if len(self.info)==1:
            input_smiles = self.info[0]["seed_smiles"]
        else:
            if self._is_not_saturated():
                input_smiles = self.info[-1]["input_smiles"]
            else:
                results = self.results[self.results["round"]==self.info[-1]["round_number"]]
                input_smiles = results['sampled_smiles'].sample().iloc[0]
        return input_smiles