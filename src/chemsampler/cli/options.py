"""Options shared by several commands."""

import sys

import click

#: The project a command works on. Everything else about where things are is in
#: its `chemsampler.json`, written by `chemsampler init`.
project_option = click.option(
    "-p",
    "--path",
    "path",
    default=None,
    type=click.Path(exists=True, file_okay=False),
    help="The project (any folder inside it). Default: the one you are in.",
)


def can_ask() -> bool:
    """Whether a person can answer a question (not a script or a pipe)."""
    return sys.stdin.isatty()


def split_ids(ids) -> list[str]:
    """Ids given with spaces, commas or both, each once, in order."""
    pieces = [piece.strip() for arg in ids for piece in arg.split(",")]
    return list(dict.fromkeys(piece for piece in pieces if piece))


def confirm(question: str, refusal: str) -> bool:
    """Ask a yes/no question, or refuse (ValueError) when nobody can answer."""
    if not can_ask():
        raise ValueError(refusal)
    return click.confirm(question, default=False)
