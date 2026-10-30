"""Options shared by several commands."""

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
