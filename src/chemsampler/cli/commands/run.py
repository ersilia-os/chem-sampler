import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("run")
@click.option(
    "--config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="JSON file with the paths and the similarity, generators and annotators "
    "lists. See the README for its keys. Any input not given by a flag, the "
    "config file or --path-to-models is looked for in the working directory.",
)
@click.option(
    "--path-to-models",
    default=None,
    type=click.Path(file_okay=False),
    help="Folder with similarity.csv, generators.csv and annotators.csv. "
    "Overrides the config file.",
)
@click.option(
    "--path-to-model-repositories",
    default=None,
    type=click.Path(file_okay=False),
    help="Folder with one sub-folder per model. Required by the run_sh backend. "
    "Overrides the config file.",
)
@click.option(
    "--similarity",
    "similarity_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with a 'similarity' column (one model). Overrides the config file "
    "and --path-to-models.",
)
@click.option(
    "--generators",
    "generators_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with a 'generators' column. Overrides the config file and "
    "--path-to-models.",
)
@click.option(
    "--annotators",
    "annotators_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with columns source, model, cutoff, direction, priority[, weight, "
    "column]. Overrides the config file and --path-to-models.",
)
@click.option(
    "--mode",
    type=click.Choice(["joint", "sequential", "weighted", "incremental"]),
    required=True,
    help="How the primary annotators combine into a round's winner.",
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
@click.option(
    "--n-steps",
    default=None,
    type=int,
    help="Number of interpolation levels per annotator. Required for "
    "--mode incremental; unused (and rejected) otherwise.",
)
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
    config: str | None,
    path_to_models: str | None,
    path_to_model_repositories: str | None,
    similarity_file: str | None,
    generators_file: str | None,
    annotators_file: str | None,
    mode: str,
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    n_steps: int | None,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: str,
    backend: str,
    output_dir: str,
    verbose: bool,
) -> None:
    """Generate and score candidate molecules across rounds, saving the result."""
    import time

    from ...config import load_inputs
    from ...optimize import hill_climb, write_results
    from ...utils.logging import logger

    if not verbose:
        logger.set_quiet_mode(True)

    start_time = time.time()
    try:
        inputs = load_inputs(
            config,
            path_to_models,
            path_to_model_repositories,
            backend=backend,
            similarity=similarity_file,
            generators=generators_file,
            annotators=annotators_file,
        )
        summary, candidates_by_round = hill_climb(
            inputs.generators,
            inputs.annotators,
            mode=mode,
            similarity=inputs.similarity,
            seed_smiles=seed_smiles,
            original_seed_smiles=original_seed_smiles,
            n_rounds=n_rounds,
            n_steps=n_steps,
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
