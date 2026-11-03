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
    help="With --model: ersilia (a Hub model) or rdkit.",
)
@click.option(
    "--model",
    "model_id",
    default=None,
    help="Score with this one model instead (a Hub id, or mw / qed / tanimoto).",
)
@click.option("--smiles", required=True, help="SMILES of the molecule to score.")
def annotate_cmd(
    path: str | None,
    source: str | None,
    model_id: str | None,
    smiles: str,
) -> None:
    """Score one molecule against the project's annotators, or one model."""
    single = source is not None or model_id is not None
    if single and (source is None or model_id is None):
        click.secho("--source and --model go together.", fg="red")
        sys.exit(1)

    from ...config import load_annotator, load_annotators, load_seeds
    from ...models.annotator import NEEDS_SEED, TanimotoAnnotator
    from ...optimize import annotate
    from ...project import PROJECT_FILE, locate_project

    # Tanimoto annotators not counted: no seed to compare with, or several
    # (one value per seed, shown apart)
    apart, per_seed = [], {}
    try:
        if single:
            # without --path, load_annotator looks for a project from here itself
            marker = None if path is None else locate_project(path) / PROJECT_FILE
            config = None if marker is None else str(marker)
            annotator = load_annotator(source, model_id, config=config)
            values = {model_id: annotator}
            if isinstance(annotator, TanimotoAnnotator):
                seeds = load_seeds(config)
                if not seeds:
                    raise ValueError(NEEDS_SEED)
                values = _per_seed(model_id, annotator, seeds)
            values = {
                label: a.score([smiles]).get(smiles, float("nan"))
                for label, a in values.items()
            }
        else:
            marker = locate_project(path) / PROJECT_FILE
            annotators = load_annotators(str(marker))
            seeds = load_seeds(str(marker))
            if len(seeds) != 1:
                apart = [
                    s for s in annotators if isinstance(s.annotator, TanimotoAnnotator)
                ]
                annotators = [s for s in annotators if s not in apart]
                if not annotators and not seeds:
                    raise ValueError(NEEDS_SEED)
                for spec in apart:
                    for label, a in _per_seed(
                        spec.annotator_id, spec.annotator, seeds
                    ).items():
                        per_seed[label] = a.score([smiles]).get(smiles, float("nan"))
            seed = seeds[0] if len(seeds) == 1 else None
            result = annotate(smiles, annotators, seed) if annotators else None
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"smiles: {smiles}")
    if single:
        for label, value in values.items():
            click.echo(f"  {label}: {value}")
        return
    row = None if result is None else result.iloc[0]
    for spec in annotators:
        click.echo(
            f"  {spec.annotator_id} [{spec.priority}] "
            f"({_describe(spec)}): {row[spec.annotator_id]}"
        )
    for spec in apart:
        lines = [
            f"{label}: {value}"
            for label, value in per_seed.items()
            if label.startswith(f"{spec.annotator_id} ")
        ]
        why = "; ".join(lines) if lines else "needs a seed (input/smiles.csv)"
        click.secho(
            f"  {spec.annotator_id} [{spec.priority}] ({_describe(spec)}): {why}; "
            "not counted",
            fg="yellow",
        )
    if row is None:
        return
    primaries = [spec for spec in annotators if spec.priority == "primary"]
    guards = [spec for spec in annotators if spec.priority == "guard-rail"]
    satisfied = int(row["cutoffs_satisfied"])
    color = "green" if satisfied == len(primaries) else "yellow"
    click.secho(f"{satisfied} of {len(primaries)} primary cutoffs satisfied", fg=color)
    if guards:
        met = bool(row["guardrails_satisfied"])
        click.secho(
            f"guard-rails {'all satisfied' if met else 'NOT all satisfied'}",
            fg="green" if met else "yellow",
        )


def _per_seed(annotator_id: str, annotator, seeds) -> dict:
    """The annotator bound to each seed, labelled by it when there are several."""
    if len(seeds) == 1:
        return {annotator_id: annotator.against(seeds[0])}
    return {
        f"{annotator_id} vs seed {k}": annotator.against(seed)
        for k, seed in enumerate(seeds, start=1)
    }


def _describe(spec) -> str:
    """The cutoff as the user wrote it: '>= 0.5', '<= 3', or '250-450'."""
    if spec.direction == "range":
        low, high = spec.cutoff
        return f"{low:g}-{high:g}"
    return f"{'>=' if spec.direction == 'higher' else '<='} {spec.cutoff:g}"
