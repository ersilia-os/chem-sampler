import sys

import click

from ..create_cli import chemsampler_cli


@chemsampler_cli.command("annotate")
@click.option(
    "--config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="JSON file with the annotators. Only its annotators and "
    "path_to_model_repositories are used. Not with --source/--model. With "
    "nothing given, ./annotators.csv in the working directory is used.",
)
@click.option(
    "--path-to-models",
    default=None,
    type=click.Path(file_okay=False),
    help="Folder with annotators.csv. Not with --source/--model.",
)
@click.option(
    "--annotators",
    "annotators_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with the annotators. Overrides the config file and "
    "--path-to-models. Not with --source/--model.",
)
@click.option(
    "--source",
    type=click.Choice(["ersilia", "rdkit"]),
    default=None,
    help="With --model: score against one annotator only, no cutoff.",
)
@click.option(
    "--model",
    "model_id",
    default=None,
    help="With --source: an Ersilia Hub model id, or an in-house name (mw, qed). "
    "Prints its raw value only.",
)
@click.option(
    "--path-to-model-repositories",
    default=None,
    type=click.Path(file_okay=False),
    help="Folder with one sub-folder per model. Required by the run_sh backend "
    "for Ersilia models.",
)
@click.option("--smiles", required=True, help="SMILES of the molecule to score.")
@click.option(
    "--backend",
    type=click.Choice(["ersilia", "run_sh"]),
    default="run_sh",
    show_default=True,
)
def annotate_cmd(
    config: str | None,
    path_to_models: str | None,
    annotators_file: str | None,
    source: str | None,
    model_id: str | None,
    path_to_model_repositories: str | None,
    smiles: str,
    backend: str,
) -> None:
    """Score a single molecule against a set of annotators, or a single model."""
    single = source is not None or model_id is not None
    a_set = (
        config is not None or path_to_models is not None or annotators_file is not None
    )
    if single and a_set:
        click.secho(
            "Use either --source/--model or --annotators/--config/--path-to-models, "
            "not both.",
            fg="red",
        )
        sys.exit(1)
    if single and (source is None or model_id is None):
        click.secho("--source and --model go together.", fg="red")
        sys.exit(1)

    from ...config import load_annotator, load_annotators
    from ...optimize import annotate

    try:
        if single:
            annotator = load_annotator(
                source, model_id, path_to_model_repositories, backend=backend
            )
            value = annotator.score([smiles]).get(smiles, float("nan"))
        else:
            annotators = load_annotators(
                config,
                path_to_models,
                path_to_model_repositories,
                backend=backend,
                annotators=annotators_file,
            )
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
