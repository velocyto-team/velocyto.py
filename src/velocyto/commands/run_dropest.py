from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from ._run import _run
from .common import init_logger, logicType, loomdtype

app = typer.Typer(name="velocyto-dropest", help="Run velocity analysis on DropEst data")


@app.callback(invoke_without_command=True)
@app.command(help="Runs the velocity analysis on DropEst preprocessed data")
def run_dropest(
    bamfile: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    gtffile: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    bcfile: Optional[Path] = typer.Option(
        None,
        "-b",
        "--bcfile",
        help="Valid barcodes file, to filter the bam. If --bcfile is not specified the file will be searched in the default position outputted by ``velocyto tools dropest_bc_correct``. Otherwise an error will be thrown",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    logic: logicType = typer.Option(
        logicType.Permissive10X, "-l", "--logic", help="The logic to use for the filtering"
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
        "--sampleid",
        "-e",
        help="The sample name that will be used as a the filename of the output.",
        exists=False,
    ),
    repmask: Optional[Path] = typer.Option(
        None,
        "-m",
        "--repmask",
        help=".gtf file containing intervals to mask ",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    samtools_threads: int = typer.Option(
        16,
        "-@",
        "--samtools-threads",
        help="The number of threads to use to sort the bam by cellID file using samtools",
    ),
    samtools_memory: str = typer.Option(
        "4G",
        "--samtools-memory",
        help="The amount of memory for samtools for each sorting thread. Accepts the same forms as samtools, so use # with K/M/G suffix",
    ),
    dtype: loomdtype = typer.Option(
        loomdtype.uint32,
        "-t",
        "--dtype",
        help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation",
    ),
    dump: str = typer.Option(
        "0",
        "-d",
        "--dump",
        help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed",
    ),
    verbose: int = typer.Option(
        0,
        "-v",
        "--verbose",
        help="Set the vebosity level: -v (only warnings) -vv (warnings and info) -vvv (warnings, info and debug)",
        count=True,
    ),
    additional_ca: typer.Context = typer.Option(...),
) -> None:
    """Runs the velocity analysis on DropEst preprocessed data

    BAMFILE  bam files to be analyzed

    GTFFILE genome annotation file
    """

    init_logger(verbose)

    additional_ca = {additional_ca[(i * 2)]: additional_ca[(i * 2) + 1] for i in range(len(additional_ca) // 2)}

    if bcfile is None:
        parentpath = bamfile.parent
        bamfilename = bamfile.name
        bcfile = parentpath.joinpath(f"barcodes_{bamfilename.split('_')[0]}.tsv")
        logger.info(f"Attempting to find automatically the valid barcode list file {bcfile}")
        if bcfile.exists():
            logger.info(f"{bcfile} found ")
        else:
            logger.info(f"{bcfile} not found!")
            logger.error(
                "In ``run_dropest`` specifying the ``--bcfile/-b`` is required. Use ``run`` if you want to adventure in a more custom usage."
            )
            logger.info("Exit without doing nothing")
            return

    if "correct" not in bamfile:
        logger.warning(
            "The file you are using does not start with the prefix ``correct_`` so it might not be the output of ``velocyto tools dropest_bc_correct``."
        )
        logger.info(
            "The program will run despite the warning but be aware of the possible consequences of not correcting the barcodes"
        )
    return _run(
        bamfile=(bamfile,),
        gtffile=gtffile,
        bcfile=bcfile,
        outputfolder=outputfolder,
        sampleid=sampleid,
        metadatatable=None,
        repmask=repmask,
        onefilepercell=False,
        logic=logic,
        without_umi=False,
        umi_extension="chr",
        multimap=False,
        test=False,
        samtools_threads=samtools_threads,
        samtools_memory=samtools_memory,
        loom_numeric_dtype=str(dtype).split(".")[-1],
        dump=dump,
        verbose=verbose,
        additional_ca=additional_ca,
    )
