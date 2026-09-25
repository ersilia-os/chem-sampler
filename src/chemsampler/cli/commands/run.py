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
    type=click.Choice(["joint", "sequential", "weighted"]),
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
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Show detailed logs with timestamps and all model calls. Default: clean progress output.",
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
    verbose: bool,
) -> None:
    """Generate and score candidate molecules across rounds, saving the result."""
    import time

    from ...config import load_annotators, load_generators
    from ...optimize import hill_climb, write_results
    from ...utils.logging import logger

    if not verbose:
        logger.set_quiet_mode(True)

    start_time = time.time()
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
        elapsed = time.time() - start_time
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)
    except KeyboardInterrupt:
        click.secho("\nInterrupted! Saving partial results...", fg="yellow")
        write_results(summary, candidates_by_round, output_dir)
        click.secho(f"Partial results written to {output_dir}", fg="yellow")
        sys.exit(0)

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

        if not verbose:
            click.echo()
            runtime_mins = elapsed / 60
            click.secho(f"Total runtime: {runtime_mins:.1f} min", fg="cyan")

            top_n = min(5, len(summary))
            top_scores = (
                summary.nlargest(top_n, "score")
                if "score" in summary.columns
                else summary.head(top_n)
            )
            if len(top_scores) > 0:
                click.secho("Top scores:", fg="cyan")
                for idx, row in top_scores.iterrows():
                    click.secho(
                        f"  Round {int(row['round'])}: {row['score']:.2f}",
                        fg="cyan",
                    )

    click.secho(f"Results written to {output_dir}", fg="green")
