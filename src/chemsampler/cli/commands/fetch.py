import dataclasses
import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option

#: The input lists `fetch IDS` writes; annotators.csv stays the user's.
_LISTS = ("generators", "similarity")


def _can_ask() -> bool:
    """Whether a person can answer a question (not a script or a pipe)."""
    return sys.stdin.isatty()


def _split(ids) -> list[str]:
    """Ids given with spaces, commas or both, each once, in order."""
    pieces = [piece.strip() for arg in ids for piece in arg.split(",")]
    return list(dict.fromkeys(piece for piece in pieces if piece))


def _listed(path) -> list[str] | None:
    """The ids a one-column list holds (none when missing), or None if it is not one."""
    try:
        lines = [line.strip() for line in path.read_text().splitlines()]
    except FileNotFoundError:
        return []
    if not lines or lines[0] != path.stem:
        return None
    return [line for line in lines[1:] if line]


@chemsampler_cli.command(
    "fetch", short_help="Clone models, download checkpoints, build their conda envs."
)
@project_option
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Replace the project's generator and similarity lists without asking.",
)
@click.argument("ids", nargs=-1)
def fetch_cmd(path: str | None, force: bool, ids: tuple[str, ...]) -> None:
    """Fetch these models (default: what input/*.csv lists), or update them."""
    from ...config import plan_fetch
    from ...hub.catalog import sort_by_role
    from ...hub.fetch import check_prerequisites, fetch_models, summarize
    from ...models.annotator import RDKIT_ANNOTATORS
    from ...project import PROJECT_FILE, Project, locate_project

    try:
        root = locate_project(path)
        roles = None
        if ids:
            roles = sort_by_role(_split(ids), in_house=tuple(RDKIT_ANNOTATORS))
            _write_lists(Project.load(root).input, roles, force)
        plan = plan_fetch(str(root / PROJECT_FILE))
        check_prerequisites()
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    if roles is not None:
        plan = _named(plan, roles, RDKIT_ANNOTATORS)
    also = ", and the rdkit environment" if plan.rdkit else ""
    click.echo(
        f"Fetching {len(plan.model_ids)} model(s): {', '.join(plan.model_ids)}{also}"
    )
    click.echo(f"  clones:       {plan.models_dir}")
    click.echo(f"  environments: {plan.envs_dir}")
    click.echo(f"  logs:         {plan.logs_dir}")
    # clones are pulled (fast-forward only) and environments that differ from
    # their recipe rebuilt: a fetched model is a current one
    results = fetch_models(plan, update=True, rebuild=True, say=click.echo)

    click.echo()
    for line in summarize(results):
        click.echo(line)
    failed = [r.model_id for r in results if not r.ok]
    if failed:
        click.secho(f"{len(failed)} model(s) failed: {', '.join(failed)}", fg="red")
        sys.exit(1)
    differing = [r.model_id for r in results if r.env == "differs"]
    if differing:
        click.secho(
            "Done, with environments that differ from their recipe.", fg="yellow"
        )
    else:
        click.secho("All models are ready.", fg="green")


def _write_lists(input_dir, roles: dict, force: bool) -> None:
    """Set generators.csv and similarity.csv to the ids, asking before a change."""
    changes = {}
    for name in _LISTS:
        path = input_dir / f"{name}.csv"
        before = _listed(path)
        if before != roles[name]:
            changes[path] = before
    if not changes:
        return
    for path, before in changes.items():
        was = "not a list" if before is None else ", ".join(before) or "none"
        now = ", ".join(roles[path.stem]) or "none (file removed)"
        click.echo(f"{path.name}: {was} -> {now}")
    if not force:
        if not _can_ask():
            raise ValueError(
                "fetch would replace the list(s) above; run it with -f to replace "
                "them without asking"
            )
        if not click.confirm("Replace?", default=False):
            click.echo("Left the lists as they are; nothing fetched.")
            sys.exit(0)
    for path in changes:
        ids = roles[path.stem]
        if ids:
            path.write_text(f"{path.stem}\n" + "".join(f"{i}\n" for i in ids))
        else:  # a run reads a missing list as none; an empty one is an error
            path.unlink(missing_ok=True)


def _named(plan, roles: dict, in_house):
    """Fetch the models named, and say where they and annotators.csv disagree."""
    samplers = [*roles["generators"], *roles["similarity"]]
    named = [m for m in roles["annotators"] if m not in in_house]
    listed = [m for m in plan.model_ids if m not in samplers]  # annotators.csv
    for model_id in named:
        if model_id not in listed:
            click.secho(
                f"{model_id} is not in input/annotators.csv: add a row with its "
                "cutoff, direction and priority before `chemsampler run`.",
                fg="yellow",
            )
    for model_id in listed:
        if model_id not in named:
            click.secho(
                f"input/annotators.csv also lists {model_id}, which you did not "
                "name: it is not fetched; name it too, or remove its row.",
                fg="yellow",
            )
    rdkit = plan.rdkit or any(m in in_house for m in roles["annotators"])
    return dataclasses.replace(plan, model_ids=(*samplers, *named), rdkit=rdkit)
