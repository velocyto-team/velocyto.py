from pathlib import Path
from typing import Optional

import typer

from ._run import _run
from .common import UMIExtension, init_logger, logicType, loomdtype

app = typer.Typer(name="velocyto-run", help="Run velocity analysis")


@app.callback(invoke_without_command=True)
@app.command(name="run")
def run(
    bamfile: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    gtffile: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    bcfile: Optional[Path] = typer.Option(
        None,
        "-b",
        "--bcfile",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Valid barcodes file, to filter the bam. If bcfile is not specified all the cell barcodes will be included. Cell barcodes should be specified in the bcfile as the `CB` tag for each read",
    ),
    outputfolder: Optional[Path] = typer.Option(
        None,
        "-o",
        "--outputfolder",
        exists=False,
        help="Output folder. If it does not exist it will be created",
    ),
    sampleid: Optional[Path] = typer.Option(
        None,
        "-e",
        "--sampleid",
        exists=False,
        help="The sample name used to retrieve informations from metadatatable",
    ),
    metadatatable: Optional[Path] = typer.Option(
        None,
        "-s",
        "--metadatatable",
        help="Table containing metadata of the various samples (csv formatted, rows are samples and cols are entries)",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    mask: Optional[Path] = typer.Option(
        None,
        "-m",
        "--mask",
        help=".gtf file containing intervals to mask",
        resolve_path=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    onefilepercell: bool = typer.Option(
        False,
        "-c",
        "--onefilepercell",
        help="If this flag is used every bamfile passed is interpreted as an independent cell, otherwise multiple files are interpreted as batch of different cells to be analyzed together. Important: cell reads distributed over multiple bamfiles is not supported!",
        is_flag=True,
    ),
    logic: logicType = typer.Option(
        logicType.Permissive10X,
        "-l",
        "--logic",
        help="Logic to use for the filtering",
    ),
    without_umi: str = typer.Option(
        False,
        "-U",
        "--without-umi",
        help="If this flag is used the data is assumed UMI-less and reads are counted instead of molecules",
        is_flag=True,
    ),
    umi_extension: UMIExtension = typer.Option(
        UMIExtension.no,
        "-u",
        "--umi-extension",
        help="In case UMI is too short to guarantee uniqueness (without information from the ampping) set this parameter to `chr`, `Gene` or `[N]bp` If set to `chr` the mapping position (binned to 10Gb intervals) will be appended to `UB` (ideal for InDrops+dropEst). If set to `Gene` then the `GX` tag will be appended to the `UB` tag. If set to `[N]bp` the first N bases of the sequence will be used to extend `UB` (ideal for STRT).",
    ),
    multimap: bool = typer.Option(
        False,
        "-M",
        "--multimap",
        help="Consider non-unique mappings (not reccomended)",
        is_flag=True,
    ),
    samtools_threads: int = typer.Option(
        16,
        "-@",
        "--samtools-threads",
        help="Number of threads to use to sort the bam by cellID file using samtools",
    ),
    samtools_memory: str = typer.Option(
        "4G",
        "--samtools-memory",
        help="The amount of memory for samtools for each sorting thread. Accepts the same forms as samtools, so use # with K/M/G suffix",
    ),
    dtype: str = typer.Option(
        loomdtype.uint32,
        "-t",
        "--dtype",
        help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation (default run: uint32)",
    ),
    dump: str = typer.Option(
        "0",
        "-d",
        "--dump",
        help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed",
        is_flag=True,
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
    """Run velocity analysis

    BAMFILE bam file with sorted reads

    GTFFILE genome annotation file
    """

    init_logger(verbose)

    additional_ca = {additional_ca[(i * 2)]: additional_ca[(i * 2) + 1] for i in range(len(additional_ca) // 2)}

    return _run(
        bamfile=bamfile,
        gtffile=gtffile,
        bcfile=bcfile,
        outputfolder=outputfolder,
        sampleid=sampleid,
        metadatatable=metadatatable,
        repmask=mask,
        onefilepercell=onefilepercell,
        logic=logic,
        without_umi=without_umi,
        umi_extension=umi_extension,
        multimap=multimap,
        test=False,
        samtools_threads=samtools_threads,
        samtools_memory=samtools_memory,
        dump=dump,
        loom_numeric_dtype=dtype.split(".")[-1],
        verbose=verbose,
        additional_ca=additional_ca,
    )
