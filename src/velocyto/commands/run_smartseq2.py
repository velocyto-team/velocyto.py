from pathlib import Path
from typing import Optional, Annotated

import typer

from velocyto.commands._run import _run
from velocyto.commands.common import init_logger, loomdtype

app = typer.Typer(
    name="velocyto-smartseq2",
    help="Run velocity analysis on SmartSeq2 data",
    rich_markup_mode="markdown",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
@app.command()
def run_smartseq2(
    bamfiles: Annotated[
        tuple[Path],
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    gtffile: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    outputfolder: Annotated[
        Optional[Path],
        typer.Option(
            "-o",
            "--outputfolder",
            help="Output folder, if it does not exist it will be created.",
            exists=False,
        ),
    ] = None,
    sampleid: Annotated[
        Optional[Path],
        typer.Option(
            "-e",
            "--sampleid",
            help="The sample name that will be used as a the filename of the output.",
            exists=False,
        ),
    ] = None,
    repmask: Annotated[
        Optional[Path],
        typer.Option(
            "-m",
            "--repmask",
            help=".gtf file containing intervals to mask",
            resolve_path=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    dtype: Annotated[
        loomdtype,
        typer.Option(
            "-t",
            "--dtype",
            help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation",
        ),
    ] = loomdtype.uint32,
    dump: Annotated[
        str,
        typer.Option(
            "-d",
            "--dump",
            help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed",
            is_flag=True,
        ),
    ] = "0",
    verbose: Annotated[
        int,
        typer.Option(
            "-v",
            "--verbose",
            help="Set the verbosity level: -v (only warnings) -vv (warnings and info) -vvv (warnings, info and debug)",
            count=True,
        ),
    ] = 0,
    additional_ca: Annotated[tuple[str], typer.Option()] = None,
) -> None:
    """Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)

    [BAMFILES, ...] a sequence of bam files to be analyzed (e.g. use a wild-card expansion)

    GTFFILE genome annotation file
    """

    init_logger(verbose)

    additional_ca = {additional_ca[(i * 2)]: additional_ca[(i * 2) + 1] for i in range(len(additional_ca) // 2)}

    return _run(
        bam_input=(bamfiles,),
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
        samtools_memory="4G",
        dump=dump,
        loom_numeric_dtype=str(dtype).split(".")[-1],
        verbose=verbose,
        additional_ca=additional_ca,
    )
