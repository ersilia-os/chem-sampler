import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option


@chemsampler_cli.command(
    "annotate", short_help="Score one molecule against the annotators (no generation)."
)
@project_option
@click.option(
    "--source",
    type=click.Choice(["ersilia", "rdkit"]),
    default=None,
    help="With --model: score against one annotator only, no cutoff. Needs no project.",
)
@click.option(
    "--model",
    "model_id",
    default=None,
    help="With --source: an Ersilia Hub model id, or an in-house name (mw, qed). "
    "Prints its raw value only.",
)
@click.option("--smiles", required=True, help="SMILES of the molecule to score.")
@click.option(
    "--backend",
    type=click.Choice(["ersilia", "run_sh"]),
    default="run_sh",
    show_default=True,
)
def annotate_cmd(
    path: str | None,
    source: str | None,
    model_id: str | None,
    smiles: str,
    backend: str,
) -> None:
    """Score a single molecule against the project's annotators, or a single model.

    With --source/--model no project is needed: the model is looked for in the
    project's folders when there is one, otherwise in this machine's store.
    """
    single = source is not None or model_id is not None
    if single and (source is None or model_id is None):
        click.secho("--source and --model go together.", fg="red")
        sys.exit(1)

    from ...config import load_annotator, load_annotators
    from ...optimize import annotate
    from ...project import PROJECT_FILE, locate_project

    try:
        if single:
            # without --path, load_annotator looks for a project from here itself
            marker = None if path is None else locate_project(path) / PROJECT_FILE
            annotator = load_annotator(
                source,
                model_id,
                backend=backend,
                config=None if marker is None else str(marker),
            )
            value = annotator.score([smiles]).get(smiles, float("nan"))
        else:
            marker = locate_project(path) / PROJECT_FILE
            annotators = load_annotators(str(marker), backend=backend)
            result = annotate(smiles, annotators)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"smiles: {smiles}")
    if single:
        click.echo(f"  {model_id}: {value}")
        return

    row = result.iloc[0]
    for spec in annotators:
        click.echo(
            f"  {spec.annotator_id} [{spec.priority}] "
            f"({_describe(spec)}): {row[spec.annotator_id]}"
        )
    primaries = [spec for spec in annotators if spec.priority == "primary"]
    guards = [spec for spec in annotators if spec.priority == "secondary"]
    satisfied = int(row["cutoffs_satisfied"])
    color = "green" if satisfied == len(primaries) else "yellow"
    click.secho(f"{satisfied} of {len(primaries)} primary cutoffs satisfied", fg=color)
    if guards:
        met = bool(row["guardrails_satisfied"])
        click.secho(
            f"guard-rails {'all satisfied' if met else 'NOT all satisfied'}",
            fg="green" if met else "yellow",
        )


def _describe(spec) -> str:
    """The cutoff as the user wrote it: '>= 0.5', '<= 3', or '250-450'."""
    if spec.direction == "range":
        low, high = spec.cutoff
        return f"{low:g}-{high:g}"
    return f"{'>=' if spec.direction == 'higher' else '<='} {spec.cutoff:g}"
