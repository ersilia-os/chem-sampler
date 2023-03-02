from chemsampler.sampler import ChemSampler
import argparse
import pandas as pd
import csv

def read_csv(in_file):
    with open(in_file, "r") as f:
        reader = csv.reader(f, delimiter=",")
        smiles_list = list(reader)
    return smiles_list[0]
    
def write_csv(out_file, smiles_list, sampled_smiles):
    dict_list = []
    for i in range(len(smiles_list)):
        smiles_dict = {}
        smiles_dict['SMILES']= smiles_list[i] 
        smiles_dict['sampled_SMILES']= sampled_smiles[i]
        dict_list.append(smiles_dict)
    field_names = ['SMILES', 'sampled_SMILES']
    with open(out_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = field_names)
        writer.writeheader()
        writer.writerows(dict_list)

def get_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("INPUT_FILE", type=str, help="Input file.")
    parser.add_argument("OUTPUT_FILE", type=str, help="Output file.")
    parser.add_argument("--SAMPLER", type=str, help="Sampler to use.")
    parser.add_argument("--NUM_SAMPLES", type=int, help="Number of samples to generate.")
    parser.add_argument("--SIM_UB", type=float, help="Upper bound of similarity.")
    parser.add_argument("--SIM_LB", type=float, help="Lower bound of similarity.")
    parser.add_argument("--DISTRIBUTION", type=str, help="Distribution of similarity.")

    return parser


def run_from_args(args: argparse.Namespace) -> None:
    sampler = ChemSampler()
    smiles_list = read_csv(args.INPUT_FILE)
    sampled_smiles = sampler.sample(
        smiles_list=smiles_list,
        Sampler=args.SAMPLER,
        num_samples=args.NUM_SAMPLES
    )
    write_csv(args.OUTPUT_FILE, smiles_list, sampled_smiles)


def main() -> None:

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main() 