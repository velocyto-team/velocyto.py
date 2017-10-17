import click
import logging
import sys
from .extract_repeats import extract_repeats
from .extract_intervals import extract_intervals
from .run import run
from .run10x import run10x
from .multi10x import multi10x


@click.group(context_settings=dict(max_content_width=300, terminal_width=300))
def cli() -> None:
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    return

cli.add_command(extract_repeats)
cli.add_command(extract_intervals)
cli.add_command(run)
cli.add_command(run10x)
cli.add_command(multi10x)
