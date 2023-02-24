import typer
from loguru import logger
from rich.console import Console
from rich.traceback import install

from . import __version__, dropest_bc_correct, run, run10x, run_dropest, run_smartseq2

install(show_locals=True, width=300, extra_lines=6, word_wrap=True)

logger.remove()

console = Console()


def version_callback(value: bool) -> None:
    """Prints the version of the package."""
    if value:
        console.print(f"[yellow]fcsjanitor[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


velocyto = typer.Typer(
    name="velocyto",
    help=("Process scRNA-seq data for RNA velocity analysis"),
    add_completion=False,
    rich_markup_mode="markdown",
)

velocyto.add_typer(run10x.app, name="run10x")
velocyto.add_typer(run.app, name="run")
velocyto.add_typer(dropest_bc_correct.app, name="dropest_bc_correct")

velocyto.add_typer(run_smartseq2.app, name="runsmartseq2")
velocyto.add_typer(run_dropest.app, name="rundropest")
