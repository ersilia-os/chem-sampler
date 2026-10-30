import os
import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option


@chemsampler_cli.command(
    "run", short_help="Generate and score candidates around a seed, over rounds."
)
@project_option
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
    "--verbose",
    is_flag=True,
    default=False,
    help="Show detailed logs with timestamps and all model calls. Default: clean progress output.",
)
def run_cmd(
    path: str | None,
    mode: str,
    seed_smiles: str | None,
    original_seed_smiles: str | None,
    n_rounds: int,
    n_steps: int | None,
    tolerance: float,
    tanimoto_cutoff: float | None,
    tanimoto_direction: str,
    backend: str,
    verbose: bool,
) -> None:
    """Generate and score candidate molecules across rounds, saving the result.

    The models and folders are the project's; each run writes a new dated folder
    under its rounds/.
    """
    import json
    import time

    import pandas as pd

    from ...config import load_inputs
    from ...optimize import hill_climb, write_results
    from ...project import PROJECT_FILE, locate_project
    from ...utils.logging import logger

    if not verbose:
        logger.set_quiet_mode(True)

    start_time = time.time()
    summary, candidates_by_round = pd.DataFrame(), {}
    output_dir = None
    try:
        root = locate_project(path)
        inputs = load_inputs(str(root / PROJECT_FILE), backend=backend)
        output_dir = str(inputs.new_round_dir())
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
        if output_dir is None:
            click.secho("\nInterrupted before the first round.", fg="yellow")
            sys.exit(0)
        click.secho("\nInterrupted! Saving partial results...", fg="yellow")
        write_results(summary, candidates_by_round, output_dir)
        click.secho(f"Partial results written to {output_dir}", fg="yellow")
        sys.exit(0)

    write_results(summary, candidates_by_round, output_dir)
    with open(os.path.join(output_dir, "inputs.json"), "w") as f:
        json.dump(inputs.describe(), f, indent=2)
        f.write("\n")

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
