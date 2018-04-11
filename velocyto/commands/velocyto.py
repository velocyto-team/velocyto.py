import click
import logging
import sys
from typing import Any
from collections import OrderedDict
from .run import run
from .run10x import run10x
from .run_smartseq2 import run_smartseq2
from .run_dropest import run_dropest
from .dropest_bc_correct import dropest_bc_correct
import velocyto._version


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


@click.version_option(version=velocyto._version.__version__)
@click.group(cls=NaturalOrderGroup, commands=OrderedDict(), context_settings=dict(max_content_width=300, terminal_width=300))
def cli() -> None:
    # NOTE: here is a good place to add a comand to control verbosity/logging-level
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    return


@click.group(cls=NaturalOrderGroup, commands=OrderedDict(), context_settings=dict(max_content_width=300, terminal_width=300))
def tools() -> None:
    """helper tools for velocyto
    """
    return

tools.add_command(dropest_bc_correct)
cli.add_command(run)
cli.add_command(run10x)
cli.add_command(run_dropest)
cli.add_command(run_smartseq2)
cli.add_command(tools)
