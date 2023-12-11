import os
import json
import pandas as pd

from chemsampler.master.master_sampler import MasterSampler
from chemsampler.utils.config import ConfigRun
from chemsampler.rules.input import InputSelector

root = os.path.dirname(os.path.abspath(__file__))
config_file = os.path.abspath(os.path.join(root, "default", "params.json"))

rc = ConfigRun(config_file)
rc.create_output_files()
params = rc.read_config_file()
rounds_info = rc.load_info_file()

inp = InputSelector(config_file)

for round in range(1, params["max_rounds"] + 1):
    input_smiles = inp.choose_input()
    ms = MasterSampler(sampler_ids=params["samplers"], descriptor_ids=params["descriptors"], unit_timeout_sec =params["time_budget_sec"])
    df, sampler_info = ms.run(seed_smiles=params["seed_smiles"], input_smiles = input_smiles, keep_smiles=params["keep_smiles"], avoid_smiles=params["avoid_smiles"])
    df["round"] = round
    df_all = rc.load_results()
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
    rc.save_results(df_updated)

    # Check if the threshold is reached
    if len(df_updated) >= params["num_samples"]:
        print(f"Threshold reached. Stopping iterations.")
        break

rc.add_info_data(rounds_info)

