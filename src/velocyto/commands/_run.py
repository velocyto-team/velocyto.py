import gzip
import itertools
import multiprocessing
import subprocess
import sys
from distutils.spawn import find_executable
from pathlib import Path
from typing import Any
from types import MappingProxyType

import joblib
import loompy
import numpy as np
import pandas as pd
import scipy as sp
from loguru import logger

from .. import version
from ..constants import BAM_COMPRESSION
from ..counter import ExInCounter
from ..metadata import MetadataCollection
from .common import choose_dtype, choose_logic, id_generator, logicType


def _run(
    *,
    bamfile: tuple[Path],
    gtffile: Path,
    bcfile: Path,
    outputfolder: Path,
    sampleid: str,
    metadatatable: str,
    repmask: str,
    onefilepercell: bool,
    logic: str,
    without_umi: str,
    umi_extension: str,
    multimap: bool,
    test: bool,
    samtools_threads: int,
    samtools_memory: str,
    loom_numeric_dtype: str,
    dump: str,
    verbose: int,
    bughunting: bool = False,
    additional_ca: dict = MappingProxyType({}),
    **kwargs,
) -> None:
    """Runs the velocity analysis outputing a loom file

    BAMFILE or [BAMFILES] one or several bam files with position-sorted

    GTFFILE annotation file

    NOTE: it is keyword only argument function
    """

    ########################
    #    Resolve Inputs    #
    ########################
    samtools = find_executable("samtools")
    if samtools is None:
        logger.error("samtools was not found")
        raise FileNotFoundError("Samtools was not found. Make sure that it is both installed and on the system path")

    loom_numeric_dtype = choose_dtype(loom_numeric_dtype)

    if isinstance(bamfile, tuple) and len(bamfile) > 1 and bamfile.suffix in [".bam", ".sam"]:
        multi = True
    elif isinstance(bamfile, tuple) and len(bamfile) == 1:
        multi = False
    else:
        raise IOError(f"Something went wrong in the argument parsing. You passed as bamfile: {bamfile}")

    if onefilepercell and multi:
        if bcfile is not None:
            raise ValueError(
                "Inputs incompatibility. --bcfile/-b option was used together with --onefilepercell/-c option."
            )
        logger.warning("Each bam file will be interpreted as a DIFFERENT cell")
    elif not onefilepercell and multi:
        logger.warning(
            "Several input files but --onefilepercell is False. Each bam file will be interpreted as containing a SET of cells!!!"
        )

    if sampleid is None:
        assert (
            metadatatable is None
        ), "--metadatatable was specified but cannot fetch sample metadata without valid sampleid"
        if multi:
            logger.warning(
                "When using mutliple files you may want to use --sampleid option to specify the name of the output file"
            )
        if multi and not onefilepercell:
            full_name = "_".join([bamfile[i].name.split(".")[0] for i in range(len(bamfile))])
            if len(full_name) > 50:
                sampleid = f'multi_input_{bamfile[0].name.split(".")[0]}_{id_generator(5)}'
            else:
                sampleid = f"multi_input_{full_name}_and_others_{id_generator(5)}"
        elif multi:
            sampleid = f'onefilepercell_{bamfile[0].name.split(".")[0]}_and_others_{id_generator(5)}'
        else:
            sampleid = f'{bamfile[0].name.split(".")[0]}_{id_generator(5)}'
        logger.info(
            f"No SAMPLEID specified, the sample will be called {sampleid} (last 5 digits are a random-id to avoid overwriting some other file by mistake)"
        )

    # Create an output folder inside the cell ranger output folder
    if outputfolder is None:
        outputfolder = bamfile[0].parent.joinpath("velocyto")
        logger.info(f"No OUTPUTFOLDER specified, find output files inside {outputfolder}")
    if not outputfolder.exists():
        outputfolder.mkdir()

    if logic not in logicType:
        raise ValueError(f"{logic} is not a valid logic. Choose one among {', '.join([_.value for _ in logicType])}")
    logger.debug(f"Using logic: {logic}")
    logic = choose_logic(logic)
    logic_obj = logic()

    if bcfile is None:
        logger.debug("Cell barcodes will be determined while reading the .bam file")
        valid_bcset = None
    else:
        # Get valid cell barcodes
        valid_bcs_list = (
            (gzip.open(bcfile).read().decode() if bcfile.suffix == ".gz" else open(bcfile).read()).rstrip().split()
        )
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])  # with sample id and with -1
        if len({bc.split("-")[0] for bc in valid_bcs_list}) == 1:
            gem_grp = f"-{valid_bcs_list[0].split('-')[-1]}"
        else:
            gem_grp = "x"
        valid_bcset = {bc.split("-")[0] for bc in valid_bcs_list}  # without -1
        logger.info(f"Read {len(valid_bcs_list)} cell barcodes from {bcfile}")
        logger.debug(f"Example of barcode: {valid_bcs_list[0].split('-')[0]} and cell_id: {valid_cellid_list[0]}")

    # Get metadata from sample sheet
    if metadatatable:
        try:
            sample_metadata = MetadataCollection(metadatatable)
            sample = sample_metadata.where("SampleID", sampleid)
            if len(sample) == 0:
                logger.error(f"Sample ID {sampleid} not found in sample sheet")
                # schema = []  # type: list
                sample = {}
            elif len(sample) > 1:
                logger.error(f"Sample ID {sampleid} has multiple lines in sample sheet")
                sys.exit(1)
            else:
                # schema = sample[0].types
                sample = sample[0].dict
            logger.debug(f"Collecting column attributes from {metadatatable}")
        except (NameError, TypeError):
            logger.warn("SAMPLEFILE was not specified. add -s SAMPLEFILE to add metadata.")
            sample = {}
    else:
        sample = {}

    ########################
    #     Start Analysis   #
    ########################

    # Initialize Exon-Intron Counter with the logic and valid barcodes (need to do it now to peek)
    if without_umi:
        if umi_extension != "no":
            logger.warning("--umi-extension was specified but incompatible with --without-umi, it will be ignored!")
        umi_extension = "without_umi"
    exincounter = ExInCounter(
        sampleid=sampleid,
        logic=logic,
        valid_bcset=valid_bcset,
        umi_extension=umi_extension,
        onefilepercell=onefilepercell,
        dump_option=dump,
        outputfolder=outputfolder,
    )

    # Heuristic to chose the memory/cpu effort
    # try:
    #     mb_available = int(subprocess.check_output("grep MemAvailable /proc/meminfo".split()).split()[1]) / 1000
    # except subprocess.CalledProcessError:
    #     logger.warning(
    #         "Your system does not support calling `grep MemAvailable /proc/meminfo` so the memory effort for the samtools command could not be chosen appropriately. 32Gb will be assumed"
    #     )
    #     mb_available = 32000  # 64Gb
    threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
    # No, I don't want anything overriding how much memory I've told samtools to use
    # mb_to_use = int(min(samtools_memory, mb_available / (len(bamfile) * threads_to_use)))
    compression = BAM_COMPRESSION

    # TODO: if we can, should check to see if the bamfile is already sorted.
    if multi and onefilepercell:
        bamfile_cellsorted = list(bamfile)
    elif onefilepercell:
        bamfile_cellsorted = [bamfile[0]]
    else:
        bamfile_cellsorted = [f"{bmf.parent.joinpath(f'cellsorted_{bmf.name}')}" for bmf in bamfile]

    if test:  # NOTE: Remove this after finishing testing, the only purpuso was to save 15min in the debugging process
        logger.warning("This place is for developer only!")

        if Path("exincounter_dump.pickle").exists():
            logger.debug("exincounter_dump.pickle is being loaded")
            exincounter = joblib.load(open("exincounter_dump.pickle", "rb"))
        else:
            logger.debug("exincounter_dump.pickle was not found")
            logger.debug("Dumping exincounter_dump.pickle BEFORE markup")
            joblib.dump(exincounter, open("exincounter_dump.pickle", "wb"))
            exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)
        check_end_process = False
    else:
        # I need to peek into the bam file to know wich cell barcode flag should be used
        if onefilepercell and without_umi:
            tagname = "NOTAG"
        elif onefilepercell:
            logger.debug("The multi input option ")
            tagname = "NOTAG"
            exincounter.peek_umi_only(bamfile[0])
        else:
            exincounter.peek(bamfile[0])
            tagname = exincounter.cellbarcode_str

        sorting_process: dict[int, Any] = {}
        for ni, bmf_cellsorted in enumerate(bamfile_cellsorted):
            # Start a subprocess that sorts the bam file
            command = f"samtools sort -l {compression} -m {samtools_memory} -t {tagname} -O BAM -@ {threads_to_use} -o {bmf_cellsorted} {bamfile[ni]}"
            if Path(bmf_cellsorted).exists():
                # This should skip sorting in smartseq2
                logger.warning(
                    f"The file {bmf_cellsorted} already exists. The sorting step will be skipped and the existing file will be used."
                )
                check_end_process = False
            else:
                sorting_process[ni] = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
                logger.info(f"Starting the sorting process of {bamfile[ni]} the output will be at: {bmf_cellsorted}")
                logger.info(f"Command being run is: {command}")
                logger.info("While the bam sorting happens do other things...")
                check_end_process = True

        # Load annotations
        logger.info(f"Load the annotation from {gtffile}")
        annotations_by_chrm_strand = exincounter.read_transcriptmodels(gtffile)
        chrs = [v for k, v in annotations_by_chrm_strand.items()]
        tms = [itertools.chain.from_iterable((v.values() for v in chrs))]
        ivls = [(itertools.chain.from_iterable(tms))]
        logger.debug(f"Generated {len(ivls)} features corresponding to {len(tms)} transcript models from {gtffile}")
        del chrs, tms, ivls

        # Load annotations
        if repmask is not None:
            logger.info(f"Load the repeat masking annotation from {repmask}")
            # mask_ivls_by_chromstrand = exincounter.read_repeats(repmask)

        # Go through the bam files a first time to markup introns
        logger.info(f"Scan {' '.join((str(_) for _ in bamfile))} to validate intron intervals")
        if bughunting and (pickled_exincounter := Path("exincounter_dump.pickle")).exists():
            with pickled_exincounter.open("rb") as exp:
                exincounter = joblib.load(exp)
        elif bughunting:
            with pickled_exincounter.open("wb") as exp:
                joblib.dump(exincounter, exp)
        else:
            exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)


    # Wait for child process to terminate
    if check_end_process:
        logger.info("Now just waiting that the bam sorting process terminates")
        for k in sorting_process.keys():
            returncode = sorting_process[k].wait()
            if returncode == 0:
                logger.info(f"bam file #{k} has been sorted")
            else:
                raise MemoryError(
                    f"bam file #{k} could not be sorted by cells.\n\
                This is probably related to an old version of samtools, please install samtools >= 1.6.\
                In alternative this could be a memory error, try to set the --samtools_memory option to a value compatible with your system. \
                Otherwise sort manually by samtools ``sort -l [compression] -m [mb_to_use]M -t [tagname] -O BAM -@ [threads_to_use] -o cellsorted_[bamfile] [bamfile]``"
                )

    # Do the actual counting
    logger.debug("Start molecule counting!")
    results = exincounter.count(
        bamfile_cellsorted, multimap=multimap
    )  # NOTE: we would avoid some millions of if statements evaluations if we write two function count and count_with output
    dict_list_arrays, cell_bcs_order = results

    ########################
    #         Output       #
    ########################

    # Prepare the loom file output
    if not exincounter.filter_mode:
        valid_bcset = exincounter.valid_bcset  # without -1
        valid_bcs_list = list(valid_bcset)  # without -1
        gem_grp = ""
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])  # with sampleid and with -1
        logger.debug(f"Example of barcode: {valid_bcs_list[0]} and cell_id: {valid_cellid_list[0]}")

    # If this data is from a 10X run, it is possible that additional_ca has data from cells that are present in the TSNE and Cluster files
    # but not present in the count matrix.   these are not removed prior to attempting to create the loom file, this whole process will fail
    if len(cell_bcs_order) < len(additional_ca["_X"]) and kwargs["is_10X"] is True:
        umap_file = list(
            Path(kwargs["samplefolder"], "outs", "per_sample_outs").rglob(
                "analysis/umap/gene_expression_2_components/projection.csv",
            )
        )[0]
        if umap_file.exists():
            umap_pd = pd.read_csv(umap_file)
            umap_pd["Barcode"] = [_[:-2] for _ in umap_pd["Barcode"]]
            umap_pd = umap_pd[umap_pd["Barcode"].isin(cell_bcs_order)]
            additional_ca["_X"] = umap_pd["UMAP-1"].to_numpy()
            additional_ca["_Y"] = umap_pd["UMAP-1"].to_numpy()

        clusters_file = list(
            Path(kwargs["samplefolder"], "outs", "per_sample_outs").rglob(
                "analysis/clustering/gene_expression_graphclust/clusters.csv",
            )
        )[0]
        if clusters_file.exists():
            # labels = np.loadtxt(clusters_file, usecols=(1,), delimiter=",", skiprows=1)
            clusters_pd = pd.read_csv(clusters_file)
            clusters_pd["Barcode"] = [_[:-2] for _ in clusters_pd["Barcode"]]
            clusters_pd = clusters_pd[clusters_pd["Barcode"].isin(cell_bcs_order)]
            additional_ca["Clusters"] = clusters_pd["Cluster"].to_numpy(dtype="int16") - 1

    ca = {
        "CellID": np.array(
            [f"{sampleid}:{v_bc}{gem_grp}" for v_bc in cell_bcs_order]
        )
    } | additional_ca
    for key, value in sample.items():
        ca[key] = np.full(len(cell_bcs_order), value)

    # Save to loom file
    outfile = outputfolder.joinpath(f"{sampleid}.loom")
    logger.debug(f"Generating output file {outfile}")

    # row attributes
    atr_table = (
        ("Gene", "genename", str),
        ("Accession", "geneid", str),
        ("Chromosome", "chrom", str),
        ("Strand", "strand", str),
        ("Start", "start", int),
        ("End", "end", int),
    )

    logger.debug("Collecting row attributes")
    ra = {}
    for name_col_attr, name_obj_attr, dtyp in atr_table:
        tmp_array = np.zeros((len(exincounter.genes),), dtype=object)  # type: np.ndarray
        for gene_id, gene_info in exincounter.genes.items():
            tmp_array[exincounter.geneid2ix[gene_id]] = getattr(gene_info, name_obj_attr)
        ra[name_col_attr] = tmp_array.astype(dtyp)

    logger.debug("Generating data table")
    layers: dict[str, np.ndarray] = {}

    for layer_name in logic_obj.layers:
        layers[layer_name] = sp.sparse.hstack(dict_list_arrays[layer_name])
        del dict_list_arrays[layer_name]

    for layer_name in logic_obj.layers:
        total: np.ndarray  # This is just a type annotation to avoid mypy complaints
        try:
            total += layers[layer_name].toarray(order="C")
        except NameError:
            total = np.array(layers[layer_name].toarray(order="C"))

    logger.debug("Writing loom file")
    # seems more likely that loompy will be some version above 2, so would be better to attempt to use the loompy2 targeting code first
    if int(loompy.__version__.split(".")[0]) >= 2:
        try:
            # tmp_layers = {"": total.astype("float32", order="C", copy=False)}
            # tmp_layers.update(
            #     {
            #         layer_name: layers[layer_name].astype(loom_numeric_dtype, copy=False)
            #         for layer_name in logic_obj.layers
            #     }
            # )
            tmp_layers = {"": total.astype("float32", order="C", copy=False)}
            # | {k: tmp_layers[k].astype(loom_numeric_dtype, copy=False) for k in tmp_layers}
            loompy.create(
                filename=str(outfile),
                layers=tmp_layers,
                row_attrs=ra,
                col_attrs=ca,
                file_attrs={
                    "velocyto.__version__": version(__name__),
                    "velocyto.logic": logic.name, # TODO: this doesn't work.  need to make string of logic type
                },
            )
            logger.debug(f"Successfully wrote to {outfile}")

        except TypeError as e:
            logger.error(e)
            sys.exit()
    elif int(loompy.__version__.split(".")[0]) < 2:
        try:
            ds = loompy.create(filename=outfile, layers=total, row_attrs=ra, col_attrs=ca)
            for layer_name in logic_obj.layers:
                ds.set_layer(name=layer_name, layers=layers[layer_name], dtype=loom_numeric_dtype)
            ds.attrs["velocyto.__version__"] = version()
            ds.attrs["velocyto.logic"] = logic
            ds.close()
        except TypeError:
            sys.exit()
    logger.debug("Terminated Succesfully!")
