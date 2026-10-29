import click


@click.group()
def chemsampler_cli() -> None:
    """chemsampler: sample and optimize chemical space around a seed molecule."""


def create_cli() -> click.Group:
    """Import every command module, registering it via its decorator's side
    effect, and return the assembled group."""
    from .commands import annotate, fetch, run, setup, update_catalog  # noqa: F401

    return chemsampler_cli
