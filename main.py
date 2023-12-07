import os
import json
import pandas as pd

from chemsampler.master_sampler import MasterSampler

root = os.path.dirname(os.path.abspath(__file__))
config_file = os.path.abspath(os.path.join(root, "default", "params.json"))

with open(config_file, 'r') as file:
    data = json.load(file)

seed_smiles = data["seed_smiles"]
keep_smiles = data["keep_smiles"]
avoid_smiles = data["avoid_smiles"]
sampler_ids = data["samplers"]
descriptor_ids = data["descriptors"]
num_samples = data["num_samples"]
max_rounds = data["max_rounds"]
time_budget_sec = data["time_budget_sec"]
output_folder = data["output_folder"]

output_folder = os.path.expanduser(output_folder)
if not os.path.exists(output_folder):
    print("path not existing")
    os.makedirs(output_folder)

output_file = os.path.join(output_folder, "chemsampler_results.csv")
info_file = os.path.join(output_folder, "rounds.json")

rounds_info = []
rounds_info += [{"seed_smiles":seed_smiles, "keep_smiles": keep_smiles, "avoid_smiles": avoid_smiles}]

for round in range(max_rounds):
    ms = MasterSampler(sampler_ids=sampler_ids, descriptor_ids = descriptor_ids, unit_timeout_sec = time_budget_sec)
    df, sampler_info = ms.run(seed_smiles=seed_smiles, input_smiles = None, keep_smiles=keep_smiles, avoid_smiles=avoid_smiles)
    df["round"] = round
    len_generated = len(df)
    if not os.path.isfile(output_file):
        df.to_csv(output_file, index=False)
    else:
        existing_data = pd.read_csv(output_file)
        # Filter out the SMILES that already exist in the masterfile
        df = df[~df["sampled_smiles"].isin(existing_data["sampled_smiles"])]
    len_not_dupl = len(df)

    round_info = {
        "round_number": round,
        "total_num_samples": len(df),
        "input_smiles": seed_smiles,
        "individual_samplers":sampler_info,
        "unique generated": len_generated - len_not_dupl
    }
    rounds_info.append(round_info)

    # Check if the threshold is reached
    if len(df) >= num_samples:
        print(f"Threshold reached. Stopping iterations.")
        break

with open(info_file, "w") as json_file:
    json.dump(rounds_info, json_file, indent=2)

