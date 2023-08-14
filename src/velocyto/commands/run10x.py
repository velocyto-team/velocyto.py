from pathlib import Path
from typing import Optional, Annotated

import typer
from loguru import logger

from velocyto.commands._run import _run
from velocyto.commands.common import init_logger, logicType, loomdtype

app = typer.Typer(
    help="Run velocity analysis on 10X Genomics data",
    rich_markup_mode="markdown",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
@app.command(
    name="run10x",
    no_args_is_help=True,
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def run10x(
    samplefolder: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=True,
            file_okay=False,
            resolve_path=True,
            readable=True,
        ),
    ],
    gtffile: Annotated[
        Path, typer.Argument(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True)
    ],
    barcode_file: Annotated[
        Optional[Path],
        typer.Option(
            "-b",
            "--barcodes",
            help="Path to the barcode file appropriate for the 10X chemistry used. Should be included with the "
            "Cellranger software under 'cellranger-X.X.X/lib/python/cellranger/barcodes/'",
            resolve_path=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    metadatatable: Annotated[
        Optional[Path],
        typer.Option(
            "-s",
            "--metadatatable",
            help="Table containing metadata of the various samples (csv fortmated rows are samples and cols are entries)",
            resolve_path=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    mask: Annotated[
        Optional[Path],
        typer.Option(
            "-m",
            "--mask",
            help=".gtf file containing intervals to mask",
            resolve_path=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    logic: Annotated[
        logicType,
        typer.Option(
            "--logic",
            "-l",
            help="The logic to use for the filtering",
        ),
    ] = logicType.Permissive10X,
    multimap: Annotated[
        bool,
        typer.Option(
            "--multimap",
            "-M",
            help="""Consider not unique mappings (not reccomended)""",
            is_flag=True,
        ),
    ] = False,
    samtools_threads: Annotated[
        int,
        typer.Option(
            "--samtools-threads",
            "-@",
            help="The number of threads to use to sort the bam by cellID file using samtools",
        ),
    ] = 16,
    samtools_memory: Annotated[
        str,
        typer.Option(
            "--samtools-memory",
            help="The amount of memory for samtools for each sorting thread. Accepts the same forms as samtools, so use # with K/M/G suffix",
        ),
    ] = "4G",
    dtype: Annotated[
        loomdtype,
        typer.Option(
            "--dtype",
            "-t",
            help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation",
        ),
    ] = loomdtype.uint16,  # why is this even an option?
    dump: Annotated[
        str,
        typer.Option(
            "-d",
            "--dump",
            help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed.",
            is_flag=True,
        ),
    ] = "0",
    verbose: Annotated[
        int,
        typer.Option(
            "--verbose",
            "-v",
            help="Set the vebosity level: -v (only warinings) -vv (warinings and info) -vvv (warinings, info and debug)",
            count=True,
        ),
    ] = 0,
    # ctx: typer.Context = typer.Option(
    #     ..., help="Extra options to pass to the job script"
    # ),
) -> None:
    """Runs the velocity analysis for a Chromium 10X Sample

    10XSAMPLEFOLDER specifies the cellranger sample folder

    GTFFILE genome annotation file
    """

    init_logger(verbose)

    # additional_ca = {ctx[(i*2)]: ctx[(i*2)+1] for i in range(len(ctx)//2)}

    # Check that the 10X analysis was run successfully
    if not samplefolder.joinpath("_log").is_file():
        logger.error(
            "This is an older version of cellranger, cannot check if the output are ready, make sure of this yourself"
        )
    elif "Pipestance completed successfully!" not in samplefolder.joinpath("_log").read_text():
        logger.error("The outputs are not ready")

    bamfile = list(samplefolder.rglob("sample_alignments.bam"))
    if not bamfile[0].exists():
        logger.error("BAM file was not found.  Are you sure you have the correct sample folder?")
        logger.exception(
            f"BAM file was not found in any subdirectories of {samplefolder}.  Are you sure you have the correct sample folder?"
        )
    elif len(bamfile) > 1:
        logger.error("Too many BAM files found. Which one?")
        logger.exception("multiple matches for `sample_alignemts.bam` were found and now I am confused. Please fix.")
    else:
        bamfile = bamfile[0]

    if not barcode_file:
        barcode_file = list(samplefolder.joinpath("outs").rglob("sample_filtered_feature_bc_matrix/barcodes.tsv.gz"))

    if not barcode_file[0].exists():
        logger.error(f"Can not locate the barcode file! Please check {barcode_file}")
    bcfile = barcode_file[0]

    outputfolder = samplefolder.joinpath("velocyto")
    sampleid = samplefolder.stem
    if outputfolder.joinpath(f"{sampleid}.loom").exists():
        raise AssertionError("The output already exist. Aborted!")

    return _run(
        bam_input=(bamfile,),
        gtffile=gtffile,
        bcfile=bcfile,
        outputfolder=outputfolder,
        sampleid=sampleid,
        metadatatable=metadatatable,
        repmask=mask,
        onefilepercell=False,
        logic=logic,
        without_umi=False,
        umi_extension="no",
        multimap=multimap,
        test=False,
        samtools_threads=samtools_threads,
        samtools_memory=samtools_memory,
        dump=dump,
        loom_numeric_dtype=str(dtype).split(".")[-1],
        verbose=verbose,
        is_10X=True,
        samplefolder=samplefolder,
    )
