from pathlib import Path
from typing import Optional

import numpy as np
import typer
from loguru import logger

from ._run import _run
from .common import init_logger, logicType, loomdtype

app = typer.Typer()  # name="velocyto-run10x", help="Run velocity analysis on 10X Genomics data")


@app.callback(invoke_without_command=True)
@app.command(
    name="run10x",
    no_args_is_help=True,
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def run10x(
    samplefolder: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    gtffile: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    metadatatable: Optional[Path] = typer.Option(
        None,
        "-s",
        "--metadatatable",
        help="Table containing metadata of the various samples (csv fortmated rows are samples and cols are entries)",
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
    logic: logicType = typer.Option(
        logicType.Permissive10X,
        "--logic",
        "-l",
        help="The logic to use for the filtering",
    ),
    multimap: bool = typer.Option(
        False,
        "--multimap",
        "-M",
        help="""Consider not unique mappings (not reccomended)""",
        is_flag=True,
    ),
    samtools_threads: int = typer.Option(
        16,
        "--samtools-threads",
        "-@",
        help="The number of threads to use to sort the bam by cellID file using samtools",
    ),
    samtools_memory: str = typer.Option(
        "4G",
        "--samtools-memory",
        help="The amount of memory for samtools for each sorting thread. Accepts the same forms as samtools, so use # with K/M/G suffix",
    ),
    dtype: loomdtype = typer.Option(
        loomdtype.uint16,
        "--dtype",
        "-t",
        help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation",
    ),  # why is this even an option?
    dump: str = typer.Option(
        "0",
        "-d",
        "--dump",
        help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed.",
        is_flag=True,
    ),
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        help="Set the vebosity level: -v (only warinings) -vv (warinings and info) -vvv (warinings, info and debug)",
        count=True,
    ),
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
        print(
            f"BAM file was not found in any subdirectories of {samplefolder}.  Are you sure you have the correct sample folder?"
        )
        exit()
    elif len(bamfile) > 1:
        logger.error("Too many BAM files found. Which one?")
        print("multiple matches for `sample_alignemts.bam` were found and now I am confused. Please fix.")
        exit()
    else:
        bamfile = bamfile[0]

    bcmatches = list(samplefolder.joinpath("outs").rglob("sample_filtered_feature_bc_matrix/barcodes.tsv.gz"))

    if not bcmatches:
        logger.error("Can not locate the barcodes.tsv file!")
    bcfile = bcmatches[0]

    outputfolder = samplefolder.joinpath("velocyto")
    sampleid = samplefolder.stem
    if outputfolder.joinpath(f"{sampleid}.loom").exists():
        raise AssertionError("The output already exist. Aborted!")

    additional_ca = {}
    try:
        umap_file = next(
            samplefolder.joinpath("outs", "per_sample_outs").rglob("umap/gene_expression_2_components/projection.csv")
        )
        if umap_file.exists():
            umap = np.loadtxt(umap_file, usecols=(1, 2), delimiter=",", skiprows=1)
            additional_ca["_X"] = umap[:, 0].astype("float32")
            additional_ca["_Y"] = umap[:, 1].astype("float32")

        clusters_file = next(
            samplefolder.joinpath("outs", "per_sample_outs").rglob("gene_expression_graphclust/clusters.csv")
        )
        if clusters_file.exists():
            labels = np.loadtxt(clusters_file, usecols=(1,), delimiter=",", skiprows=1)
            additional_ca["Clusters"] = labels.astype("int") - 1

    except Exception:
        logger.error("Some IO problem in loading cellranger umap/pca/kmeans files occurred!")

    return _run(
        bamfile=(bamfile,),
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
        additional_ca=additional_ca,
        is_10X=True,
        samplefolder=samplefolder,
    )
