import json
import os
import pandas as pd

class ConfigRun(object):
    def __init__(self, config_file):
        self.config_file = config_file
    
    def read_config_file(self):
        with open(self.config_file, "r") as json_file:
            params = json.load(json_file)
        return params
    
    def _create_output_folder(self):
        params= self.read_config_file()
        output_folder = params["output_folder"]
        output_folder = os.path.expanduser(output_folder)
        if not os.path.exists(output_folder):
            print("Output Folder Not Existing")
            os.makedirs(output_folder)
        return output_folder
    
    def _results_file(self):
        output_folder = self._create_output_folder()
        results_file = os.path.join(output_folder, "chemsampler_results.csv")
        seed_smiles = self.read_config_file()["seed_smiles"]
        df = pd.DataFrame({"round": [0], "sampled_smiles":[seed_smiles]})
        df.to_csv(results_file, index=False)    

    def _round_info_file(self):
        output_folder = self._create_output_folder()
        info_file = os.path.join(output_folder, "rounds_info.json")
        params = self.read_config_file()
        rounds_info = []
        rounds_info += [{"seed_smiles":params["seed_smiles"], "keep_smiles": params["keep_smiles"], "avoid_smiles": params["avoid_smiles"]}]
        with open(info_file, "w") as json_file:
            json.dump(rounds_info, json_file, indent=2)

    def load_results(self):
        output_folder = self._create_output_folder()
        results_file = os.path.join(output_folder, "chemsampler_results.csv")
        df = pd.read_csv(results_file)
        return df    
    
    def save_results(self, df):
        output_folder = self._create_output_folder()
        results_file = os.path.join(output_folder, "chemsampler_results.csv")
        df.to_csv(results_file, index=False)
    
    def load_info_file(self):
        output_folder = self._create_output_folder()
        info_file = os.path.join(output_folder, "rounds_info.json")
        with open(info_file, 'r') as json_file:
            rounds_info = json.load(json_file)
        return rounds_info
    
    def add_info_data(self, rounds_info):
        output_folder = self._create_output_folder()
        info_file = os.path.join(output_folder, "rounds_info.json")
        with open(info_file, "w") as json_file:
            json.dump(rounds_info, json_file, indent=2)

    def create_output_files(self):
        self._results_file()
        self._round_info_file()

