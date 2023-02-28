#!/usr/bin/env python3
import argparse
import pandas as pd
import csv
from molecule_generation import load_model_from_directory
from molecule_generation.utils.cli_utils import (
    get_model_loading_parser,
    setup_logging,
    supress_tensorflow_warnings,
)

def print_samples(model_dir: str, num_samples: int, **model_kwargs) -> None:
    with load_model_from_directory(model_dir, **model_kwargs) as model:
        samples = model.sample(num_samples)

    print("\n".join(samples))

def save_samples(model_dir: str, num_samples: int, data_file: str, **model_kwargs,) -> None:
    with load_model_from_directory(model_dir, **model_kwargs) as model:
        samples = model.sample(num_samples)
    dict = {'smiles': samples}
    df_samples = pd.DataFrame(dict)
    df_samples.to_csv(data_file, index=False)


def get_argparser() -> argparse.ArgumentParser:
    parser = get_model_loading_parser(description="Sample SMILES strings from a trained model.")
    parser.add_argument("NUM_SAMPLES", type=int, help="Number of samples to generate.")
    parser.add_argument("DATA_FILE", type=str, help="name of the output file")

    parser.add_argument(
        "--beam-size", dest="beam_size", type=int, help="Beam size to use during decoding."
    )

    return parser


def run_from_args(args: argparse.Namespace) -> None:
    model_kwargs = {key: getattr(args, key) for key in ["beam_size", "seed", "num_workers"]}

    save_samples(
        model_dir=args.MODEL_DIR,
        num_samples=args.NUM_SAMPLES,
        data_file=args.DATA_FILE,
        **{key: value for (key, value) in model_kwargs.items() if value is not None},
        
    )


def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main()
