import dataclasses
import sys

import click

from ..create_cli import chemsampler_cli
from ..options import project_option

#: The input lists `fetch IDS` writes; annotators.csv stays the user's.
_LISTS = ("generators", "similarity")

#: The id that asks for the `rdkit` environment, which mw, qed and tanimoto run in.
_RDKIT = "rdkit"


def _can_ask() -> bool:
    """Whether a person can answer a question (not a script or a pipe)."""
    return sys.stdin.isatty()


def _split(ids) -> list[str]:
    """Ids given with spaces, commas or both, each once, in order."""
    pieces = [piece.strip() for arg in ids for piece in arg.split(",")]
    return list(dict.fromkeys(piece for piece in pieces if piece))


def _pinned(name: str) -> tuple[str, str]:
    """`eos9p57@v1.0.0` -> ("eos9p57", "v1.0.0"); a bare id is `live`."""
    from ...hub.fetch import LIVE

    model_id, _, tag = name.partition("@")
    if "@" in name and not tag:
        raise ValueError(
            f"{name}: give the version after @: a release tag (e.g. "
            f"{model_id}@v1.0.0) or {LIVE}"
        )
    return model_id, tag or LIVE


def _listed(path) -> list[tuple[str, str]] | None:
    """The (model, version) rows of a model list (none when missing), or None if
    it is not one."""
    import csv

    try:
        with open(path, newline="") as f:
            rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
    except FileNotFoundError:
        return []
    if not rows or [c.strip() for c in rows[0]] != ["model", "version"]:
        return None
    return [tuple(c.strip() for c in r) for r in rows[1:]]


def _confirm(question: str, refusal: str) -> bool:
    """Ask a yes/no question, or refuse when nobody can answer."""
    if not _can_ask():
        raise ValueError(refusal)
    return click.confirm(question, default=False)


@chemsampler_cli.command(
    "fetch", short_help="Clone models, download checkpoints, build their conda envs."
)
@project_option
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Do not ask: replace the project's model lists, or fetch all they list.",
)
@click.argument("ids", nargs=-1)
def fetch_cmd(path: str | None, force: bool, ids: tuple[str, ...]) -> None:
    """Fetch these models, ID (live) or ID@TAG (default: all input/*.csv lists, asking)."""
    from ...config import plan_fetch, plan_rdkit
    from ...hub.catalog import sort_by_role
    from ...hub.fetch import check_prerequisites, fetch_models, summarize
    from ...models.annotator import RDKIT_ANNOTATORS
    from ...project import PROJECT_FILE, Project, locate_project

    try:
        pinned = dict(_pinned(name) for name in _split(ids))
        names = list(pinned)
        rdkit = _RDKIT in names
        root = locate_project(path)
        config = str(root / PROJECT_FILE)
        in_house = [n for n in names if n in RDKIT_ANNOTATORS]
        if in_house:
            raise ValueError(
                f"{', '.join(in_house)}: not models but RDKit annotators, which run "
                f"in the rdkit environment: use `chemsampler fetch {_RDKIT}`"
            )
        if names == [_RDKIT]:
            plan = plan_rdkit(config)  # the env alone: the lists stay as they are
        elif names:
            # `rdkit` stands for the annotators it runs, so a set with no other
            # annotator is still a set
            roles = sort_by_role(names, in_house=(_RDKIT,))
            roles["annotators"] = [m for m in roles["annotators"] if m != _RDKIT]
            samplers = [*roles["generators"], *roles["similarity"]]
            versions = {m: pinned[m] for m in samplers}
            _write_lists(Project.load(root).input, roles, versions, force)
            plan = _named(plan_fetch(config), roles, rdkit, pinned)
        else:
            plan = plan_fetch(config)
            plan = dataclasses.replace(plan, rdkit=bool(plan.rdkit_annotators))
            click.echo(f"input/*.csv list: {_what(plan)}")
            if not force and not _confirm(
                "Fetch these?",
                "fetch would fetch everything listed above; run it with -f to "
                "fetch without asking",
            ):
                click.echo("Nothing fetched.")
                sys.exit(0)
        check_prerequisites()
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    click.echo(f"Fetching {_what(plan)}")
    click.echo(f"  clones:       {plan.models_dir}")
    click.echo(f"  environments: {plan.envs_dir}")
    click.echo(f"  logs:         {plan.logs_dir}")
    # clones are checked out at their release and environments that differ from
    # their recipe rebuilt: a fetched model is the version the inputs name
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


def _what(plan) -> str:
    """'2 model(s): eos9taz v1.0.0, eos8nyg v1.0.0, and the rdkit environment'."""
    models = [f"{m} {plan.versions.get(m, '?')}" for m in plan.model_ids]
    env = "the rdkit environment" + (
        f" ({', '.join(plan.rdkit_annotators)})" if plan.rdkit_annotators else ""
    )
    if not models:
        return env
    also = f", and {env}" if plan.rdkit else ""
    return f"{len(models)} model(s): {', '.join(models)}{also}"


def _write_lists(input_dir, roles: dict, versions: dict, force: bool) -> None:
    """Set generators.csv and similarity.csv to the ids, asking before a change."""
    changes = {}
    for name in _LISTS:
        path = input_dir / f"{name}.csv"
        before = _listed(path)
        wanted = [(m, versions[m]) for m in roles[name]]
        if before != wanted:
            changes[path] = (before, wanted)
    if not changes:
        return
    for path, (before, wanted) in changes.items():
        was = "not a list" if before is None else _shown(before) or "none"
        now = _shown(wanted) or "none (file removed)"
        click.echo(f"{path.name}: {was} -> {now}")
    if not force and not _confirm(
        "Replace?",
        "fetch would replace the list(s) above; run it with -f to replace them "
        "without asking",
    ):
        click.echo("Left the lists as they are; nothing fetched.")
        sys.exit(0)
    for path, (_, wanted) in changes.items():
        if wanted:
            rows = "".join(f"{m},{v}\n" for m, v in wanted)
            path.write_text("model,version\n" + rows)
        else:  # a run reads a missing list as none; an empty one is an error
            path.unlink(missing_ok=True)


def _shown(pairs) -> str:
    return ", ".join(f"{m} {v}" for m, v in pairs)


def _named(plan, roles: dict, rdkit: bool, pinned: dict):
    """
    Fetch the models named, and say where they and annotators.csv disagree. An
    annotator is fetched at the version annotators.csv gives it, if it has a row.
    """
    samplers = [*roles["generators"], *roles["similarity"]]
    named = roles["annotators"]
    listed = [m for m in plan.model_ids if m not in samplers]  # annotators.csv
    versions = dict(plan.versions)
    for model_id in named:
        if model_id in listed:
            tag = pinned[model_id]
            if tag != versions[model_id]:
                click.secho(
                    f"input/annotators.csv lists {model_id} at {versions[model_id]}, "
                    f"not {tag}: fetching {versions[model_id]}; change its row to "
                    f"use {tag}.",
                    fg="yellow",
                )
            continue
        versions[model_id] = pinned[model_id]
        click.secho(
            f"{model_id} is not in input/annotators.csv: add a row with its version "
            f"({versions[model_id]}), cutoff, direction and priority before "
            "`chemsampler run`.",
            fg="yellow",
        )
    for model_id in listed:
        if model_id not in named:
            click.secho(
                f"input/annotators.csv also lists {model_id}, which you did not "
                "name: it is not fetched; name it too, or remove its row.",
                fg="yellow",
            )
    return dataclasses.replace(
        plan, model_ids=(*samplers, *named), rdkit=rdkit, versions=versions
    )
