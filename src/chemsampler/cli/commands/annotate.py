import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("annotate")
@click.option(
    "--annotators",
    "annotators_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with columns annotator_id, cutoff, direction[, column]. "
    "Mutually exclusive with --model.",
)
@click.option(
    "--model",
    "model_id",
    default=None,
    help="Single annotator id ('qed' or an Ersilia Hub model id); prints its "
    "raw value only, no cutoff. Mutually exclusive with --annotators.",
)
@click.option("--smiles", required=True, help="SMILES of the molecule to score.")
@click.option(
    "--backend",
    type=click.Choice(["ersilia", "run_sh"]),
    default="run_sh",
    show_default=True,
)
def annotate_cmd(
    annotators_path: str | None, model_id: str | None, smiles: str, backend: str
) -> None:
    """Score a single molecule against a set of annotators, or a single --model."""
    if (annotators_path is None) == (model_id is None):
        click.secho("Exactly one of --annotators or --model is required.", fg="red")
        sys.exit(1)

    from ...config import build_annotator, load_annotators
    from ...optimize import annotate

    try:
        if model_id is not None:
            annotator = build_annotator(model_id, backend=backend)
            value = annotator.score([smiles]).get(smiles, float("nan"))
        else:
            annotators = load_annotators(annotators_path, backend=backend)
            result = annotate(smiles, annotators)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"smiles: {smiles}")
    if model_id is not None:
        click.echo(f"  {model_id}: {value}")
        return

    row = result.iloc[0]
    for spec in annotators:
        click.echo(f"  {spec.annotator_id}: {row[spec.annotator_id]}")
    satisfied = int(row["cutoffs_satisfied"])
    total = len(annotators)
    color = "green" if satisfied == total else "yellow"
    click.secho(f"{satisfied} of {total} cutoffs satisfied", fg=color)
