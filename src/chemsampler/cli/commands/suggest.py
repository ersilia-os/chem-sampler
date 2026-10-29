import sys

import click

from ..create_cli import chemsampler_cli

_DIRECTIONS = ["higher", "lower"]


def _is_terminal() -> bool:
    """Whether a person is reading stdout: a table for them, plain CSV otherwise."""
    return sys.stdout.isatty()


def _percentile_grid(pool):
    """The library's values as a one-row table: min, the percentiles, max."""
    from rich import box
    from rich.table import Table

    labels = ["min", *(f"p{level:g}" for level, _ in pool.percentiles), "max"]
    values = [pool.low, *(value for _, value in pool.percentiles), pool.high]
    grid = Table(box=box.SIMPLE_HEAD, header_style="dim", pad_edge=False)
    for label in labels:
        grid.add_column(label, justify="right")
    grid.add_row(*(f"{value:.5g}" for value in values))
    return grid


def _show_table(result, pool, entry, model, version, column, direction) -> None:
    """The result for a person at a terminal."""
    from rich import box
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text

    # Text, never a bare string: Rich would read the [nH] in a SMILES as markup and
    # drop it
    console = Console(highlight=False)
    console.print(Text(entry.title, style="bold"))
    console.print(Text(f"{model} · {version} · {column}", style="dim"))
    if pool is None:
        console.print(Text("No compound has a value for this column."))
        return
    end = "highest" if direction == "higher" else "lowest"
    console.print(
        Text(
            f"{len(result)} {end} of {pool.count:,} scored compounds "
            "in the Ersilia Reference Library"
        )
    )
    console.print(_percentile_grid(pool))
    table = Table(box=box.ROUNDED, header_style="bold")
    table.add_column("#", justify="right", style="dim")
    table.add_column("score", justify="right")
    table.add_column("SMILES", overflow="fold")
    for rank, (smiles, score) in enumerate(zip(result["smiles"], result[column]), 1):
        table.add_row(str(rank), f"{score:.6g}", Text(smiles))
    console.print(table)
    console.print(Text("Save as CSV: add > seeds.csv to this command.", style="dim"))


@chemsampler_cli.command(
    "suggest", short_help="Top-N compounds from isaura precalculations (developers)."
)
@click.option(
    "--model", required=True, help="Ersilia annotation model id, e.g. eos4zfy."
)
@click.option(
    "--version",
    required=True,
    help="Stored version, e.g. v2. Versions of one model can give different results.",
)
@click.option("--column", required=True, help="Output column of the model.")
@click.option(
    "--direction",
    required=True,
    type=click.Choice(_DIRECTIONS),
    help="higher: largest values first. lower: smallest first.",
)
@click.option(
    "--top",
    default=10,
    show_default=True,
    type=click.IntRange(min=1),
    help="How many compounds.",
)
def suggest_cmd(
    model: str, version: str, column: str, direction: str, top: int
) -> None:
    """Developer tool: the compounds with the highest or lowest value of one column.

    Only for annotation models (the ones that give scores); representation and
    sampling models are refused. Reads the precalculations isaura holds for the
    Ersilia Reference Library, so no model is run. Needs `pip install
    chemsampler[isaura]` and `isaura configure --remote`. At a terminal it shows a
    table; redirected or piped it prints plain CSV (smiles, COLUMN), to keep as
    seeds for `chemsampler run`.
    """
    from ...hub.precalc import (
        PrecalcError,
        open_isaura,
        pool_stats,
        require_annotator,
        top_compounds,
    )

    terminal = _is_terminal()
    try:
        entry = require_annotator(model)
        connection, files = open_isaura(model, version)
        result = top_compounds(connection, files, column, direction, top)
        # the context costs one more pass over the column, so only a person gets it
        pool = (
            pool_stats(connection, files, column) if terminal and len(result) else None
        )
    except PrecalcError as e:
        click.secho(str(e), fg="red", err=True)
        sys.exit(1)

    if terminal:
        _show_table(result, pool, entry, model, version, column, direction)
    else:
        click.echo(result.to_csv(index=False), nl=False)
