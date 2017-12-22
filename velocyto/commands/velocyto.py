import click
import logging
import sys
from typing import Any
from collections import OrderedDict
from .run import run
from .run10x import run10x


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.

    Make sure you initialize the `self.commands` with OrderedDict instance.

    With decorator, use::

        @click.group(cls=NaturalOrderGroup, commands=OrderedDict())
    """

    def list_commands(self, ctx: Any) -> Any:
        """List command names as they are in commands dict.

        If the dict is OrderedDict, it will preserve the order commands
        were added.
        """
        return self.commands.keys()


@click.group(cls=NaturalOrderGroup, commands=OrderedDict(), context_settings=dict(max_content_width=300, terminal_width=300))
def cli() -> None:
    # NOTE: here is a good place to add a comand to control verbosity/logging-level
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    return

cli.add_command(run)
cli.add_command(run10x)
