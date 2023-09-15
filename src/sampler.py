from chemsampler.sampler import ChemSampler
import argparse
import pandas as pd
import csv


def read_csv(in_file):
    df = pd.read_csv(in_file)
    col_name = df.columns[0]
    smiles_list = df[col_name].tolist()
    return smiles_list


def write_csv(out_file, sampled_smiles, sampler):
    if sampler == "all":
        with open(out_file, "w", newline="") as file:
            writer = csv.writer(file)
            for i in range(len(sampled_smiles)):
                key = list(sampled_smiles.keys())[i]
                # writer.writerow([key]) #it will write a single list, if want to retrieve by sampler, uncomment this
                for j in range(len(sampled_smiles[key])):
                    writer.writerow([sampled_smiles[key][j]])
    else:
        with open(out_file, "w", newline="") as file:
            writer = csv.writer(file)
            for i in range(len(sampled_smiles)):
                writer.writerow([sampled_smiles[i]])


def get_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("INPUT_FILE", type=str, help="Input file.")
    parser.add_argument("OUTPUT_FILE", type=str, help="Output file.")
    parser.add_argument("--sampler", type=str, help="Sampler to use.")
    parser.add_argument(
        "--num_samples", type=int, help="Number of samples to generate."
    )
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
        sim_lb=args.sim_lb,
        distribution=args.distribution,
    )
    write_csv(args.OUTPUT_FILE, sampled_smiles, args.sampler)


def main() -> None:

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main()
