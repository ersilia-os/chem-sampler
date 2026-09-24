import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("annotate")
@click.option(
    "--annotators",
    "annotators_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with columns annotator_id, cutoff, direction[, column].",
)
@click.option("--smiles", required=True, help="SMILES of the molecule to score.")
@click.option(
    "--backend",
    type=click.Choice(["ersilia", "run_sh"]),
    default="ersilia",
    show_default=True,
)
def annotate_cmd(annotators_path: str, smiles: str, backend: str) -> None:
    """Score a single molecule against a set of annotators."""
    from ...config import load_annotators
    from ...optimize import annotate

    try:
        annotators = load_annotators(annotators_path, backend=backend)
        result = annotate(smiles, annotators)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    row = result.iloc[0]
    click.echo(f"smiles: {smiles}")
    for spec in annotators:
        click.echo(f"  {spec.annotator_id}: {row[spec.annotator_id]}")

    satisfied = int(row["cutoffs_satisfied"])
    total = len(annotators)
    color = "green" if satisfied == total else "yellow"
    click.secho(f"{satisfied} of {total} cutoffs satisfied", fg=color)
