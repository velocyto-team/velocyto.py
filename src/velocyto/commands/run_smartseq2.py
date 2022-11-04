import sys
from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from ._run import _run
from .common import loomdtype

app = typer.Typer(name="velocyto-smartseq2", help="Run velocity analysis on SmartSeq2 data")


@app.callback(invoke_without_command=True)
@app.command()
def run_smartseq2(
    bamfiles: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    gtffile: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    outputfolder: Optional[Path] = typer.Option(
        None,
        "-o",
        "--outputfolder",
        help="Output folder, if it does not exist it will be created.",
        exists=False,
    ),
    sampleid: Optional[Path] = typer.Option(
        None,
        "-e",
        "--sampleid",
        help="The sample name that will be used as a the filename of the output.",
        exists=False,
    ),
    repmask: Optional[Path] = typer.Option(
        None,
        "-m",
        "--repmask",
        help=".gtf file containing intervals to mask",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    dtype: loomdtype = typer.Option(
        loomdtype.uint32,
        "-t",
        "--dtype",
        help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation",
    ),
    dump: bool = typer.Option(
        False,
        "-d",
        "--dump",
        help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed",
        is_flag=True,
    ),
    verbose: int = typer.Option(
        0,
        "-v",
        "--verbose",
        help="Set the verbosity level: -v (only warnings) -vv (warnings and info) -vvv (warnings, info and debug)",
        count=True,
    ),
    additional_ca: typer.Context = typer.Option(...),
) -> None:
    """Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)

    [BAMFILES, ...] a sequence of bam files to be analyzed (e.g. use a wild-card expansion)

    GTFFILE genome annotation file
    """

    if verbose == 3:
        logger.add(sys.stderr, level="DEBUG")
    elif verbose == 2:
        logger.add(sys.stderr, level="INFO")
    elif verbose == 1:
        logger.add(sys.stderr, level="WARNING")
    else:
        logger.add(sys.stderr, level="ERROR")

    additional_ca = {additional_ca[(i * 2)]: additional_ca[(i * 2) + 1] for i in range(len(additional_ca) // 2)}

    return _run(
        bamfile=bamfiles,
        gtffile=gtffile,
        bcfile=None,
        outputfolder=outputfolder,
        sampleid=sampleid,
        metadatatable=None,
        repmask=repmask,
        onefilepercell=True,
        logic="SmartSeq2",
        without_umi=True,
        umi_extension="no",
        multimap=False,
        test=False,
        samtools_threads=1,
        samtools_memory=1,
        dump=dump,
        loom_numeric_dtype=dtype,
        verbose=verbose,
        additional_ca=additional_ca,
    )
