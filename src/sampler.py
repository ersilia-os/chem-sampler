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
    parser.add_argument("--sampler", type=str, help="Sampler to use.")
    parser.add_argument("--num_samples", type=int, help="Number of samples to generate.")
    parser.add_argument("--sim_ub", type=float, help="Upper bound of similarity.")
    parser.add_argument("--sim_lb", type=float, help="Lower bound of similarity.")
    parser.add_argument("--distribution", type=str, help="Distribution of similarity.")

    return parser


def run_from_args(args: argparse.Namespace) -> None:
    sampler = ChemSampler()
    smiles_list = read_csv(args.INPUT_FILE)
    sampled_smiles = sampler.sample(
        smiles_list=smiles_list,
        Sampler=args.sampler,
        num_samples=args.num_samples,
        sim_ub=args.sim_ub,
        sim_lb= args.sim_lb,
        distribution = args.distribution,
    )
    write_csv(args.OUTPUT_FILE, smiles_list, sampled_smiles)


def main() -> None:

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main() 