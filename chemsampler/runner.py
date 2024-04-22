from .master.master_sampler import MasterSampler
from .utils.config import ConfigRun
from .rules.input import InputSelector
from .utils.visualize import VisualizeMolecules


class Runner(object):
    
    def __init__(self, config_file):
        self.config_file = config_file

    def run(self):
        rc = ConfigRun(self.config_file)
        output_folder = rc._create_output_folder()
        rc.create_output_files()
        params = rc.read_config_file()

        for round in range(1, params["max_rounds"] + 1):
            rounds_info = rc.load_info_file()
            df_all = rc.load_results()
            inp = InputSelector(rounds_info, df_all, params["saturation_number"])
            input_smiles = inp.choose_input()
            ms = MasterSampler(sampler_ids=params["samplers"], descriptor_ids=params["descriptors"], unit_timeout_sec =params["time_budget_sec"])
            df, sampler_info, df_ = ms.run(seed_smiles=params["seed_smiles"], input_smiles = input_smiles, keep_smiles=params["keep_smiles"], avoid_smiles=params["avoid_smiles"])
            df["round"] = round
            len_generated = len(df)
            df = df[~df["sampled_smiles"].isin(df_all["sampled_smiles"])]
            len_not_dupl = len(df)
            df_all.set_index(["sampled_smiles", "round"], inplace=True)
            df.set_index(["sampled_smiles", "round"], inplace=True)
            df_updated = df_all.combine_first(df).fillna(0).reset_index()
            round_info = {
                "round_number": round,
                "total_num_samples": len(df_updated),
                "input_smiles": input_smiles,
                "individual_samplers":sampler_info,
                "total_generated":len_generated,
                "total_unique_generated": len_not_dupl
            }
            rounds_info.append(round_info)
            rc.add_info_data(rounds_info)
            rc.save_results(df_updated)
            rc.save_discarded(df_)
            df.reset_index(inplace=True)
            vm = VisualizeMolecules()
            vm.visualize_mols(df, output_folder, round)

            # Check if the threshold is reached
            if len(df_updated) >= params["num_samples"]:
                print(f"Threshold reached. Stopping iterations.")
                break
        # at the end of the loop, add a molecule id column
        num_digits = len(str(len(df_updated)))
        df_updated["id"] = ["chemsampler-{:0{}}".format(i, num_digits) for i in range(1, len(df_updated) + 1)]
        rc.save_results(df_updated)

    def add_properties(self):
        rc = ConfigRun(self.config_file)
        df = rc.add_calculated_properties()
        print(df.columns)
        cols = ["id", "sampled_smiles", "round", "mw", "qed", "logp"]
        other_cols = [col for col in df.columns if col not in cols]
        col_order = cols + other_cols
        df = df[col_order]
        rc.save_results(df)
