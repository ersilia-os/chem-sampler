import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("run")
@click.option(
    "--annotators",
    "annotators_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with columns annotator_id, cutoff, direction[, column]. Falls "
    "back to ./annotators.csv in the cwd.",
)
@click.option(
    "--generators",
    "generators_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with a generator_id column. Falls back to ./generators.csv in "
    "the cwd, then the 6 default Hub generators.",
)
@click.option(
    "--mode",
    type=click.Choice(["joint", "sequential"]),
    required=True,
    help="How annotators combine into a round's winner.",
)
@click.option("--seed-smiles", default=None, help="SMILES of the starting molecule.")
@click.option(
    "--original-seed-smiles",
    default=None,
    help="SMILES of the true original molecule in a manually re-seeded chain "
    "of runs. Adds a tanimoto_to_original_seed column, independent of "
    "--seed-smiles.",
)
@click.option("--n-rounds", default=5, show_default=True, type=int)
@click.option("--tolerance", default=0.0, show_default=True, type=float)
@click.option(
    "--tanimoto-cutoff",
    default=None,
    type=float,
    help="Similarity cutoff vs --seed-smiles, gated by --tanimoto-direction. "
    "Requires --seed-smiles.",
)
@click.option(
    "--tanimoto-direction",
    type=click.Choice(["higher", "lower"]),
    default="higher",
    show_default=True,
    help="'higher' keeps candidates at least this similar to the seed; "
    "'lower' pushes toward novelty (at most this similar).",
)
@click.option(
    "--backend",
    type=click.Choice(["ersilia", "run_sh"]),
    default="run_sh",
    show_default=True,
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write summary.csv and round<n>.csv into.",
)
def run_cmd(
    annotators_path: str | None,
    generators_path: str | None,
    mode: str,
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: str,
    backend: str,
    output_dir: str,
) -> None:
    """Generate and score candidate molecules across rounds, saving the result."""
    from ...config import load_annotators, load_generators
    from ...optimize import hill_climb, write_results

    try:
        annotators = load_annotators(annotators_path, backend=backend)
        generator = load_generators(generators_path, backend=backend)
        summary, candidates_by_round = hill_climb(
            generator,
            annotators,
            mode=mode,
            seed_smiles=seed_smiles,
            original_seed_smiles=original_seed_smiles,
            n_rounds=n_rounds,
            tolerance=tolerance,
            tanimoto_cutoff=tanimoto_cutoff,
            tanimoto_direction=tanimoto_direction,
        )
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    write_results(summary, candidates_by_round, output_dir)

    if summary.empty:
        click.secho("No candidates were generated; nothing to report.", fg="yellow")
    else:
        best = summary.iloc[-1]
        click.secho(
            f"Best candidate: {best['smiles']} (round {best['round']}, "
            f"score={best['score']})",
            fg="green",
        )
    click.secho(f"Results written to {output_dir}", fg="green")
