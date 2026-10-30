import click

#: The order a user needs the commands in, rather than alphabetical.
_COMMAND_ORDER = ("init", "fetch", "run", "annotate", "catalog", "suggest")

_EPILOG = """\b
Typical use:
  chemsampler init -p my-run && cd my-run
  (edit input/*.csv)
  chemsampler fetch
  chemsampler run --mode sequential --seed-smiles SMILES"""


class _WorkflowGroup(click.Group):
    """A group that lists its commands in the order of the workflow."""

    def list_commands(self, ctx: click.Context) -> list[str]:
        known = len(_COMMAND_ORDER)
        return sorted(
            super().list_commands(ctx),
            key=lambda name: (
                _COMMAND_ORDER.index(name) if name in _COMMAND_ORDER else known,
                name,
            ),
        )


@click.group(cls=_WorkflowGroup, epilog=_EPILOG)
def chemsampler_cli() -> None:
    """chemsampler: sample and optimize chemical space around a seed molecule."""


def create_cli() -> click.Group:
    """Import every command module, registering it via its decorator's side
    effect, and return the assembled group."""
    from .commands import annotate, catalog, fetch, init, run, suggest  # noqa: F401

    return chemsampler_cli
