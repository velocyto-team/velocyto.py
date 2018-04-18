import sys
import os
import glob
import re
import gzip
import array
import loompy
import numpy as np
import random
import string
import subprocess
import multiprocessing
import csv
import itertools
from collections import defaultdict
import logging
import h5py
from typing import *
import velocyto as vcy


def id_generator(size: int=6, chars: str=string.ascii_uppercase + string.digits) -> str:
    return ''.join(random.choice(chars) for _ in range(size))


def _run(*, bamfile: Tuple[str], gtffile: str,
         bcfile: str, outputfolder: str,
         sampleid: str, metadatatable: str,
         repmask: str, onefilepercell: bool, logic: str,
         without_umi: str, umi_extension: str, multimap: bool, test: bool,
         samtools_threads: int, samtools_memory: int, dump: bool, verbose: int,
         additional_ca: dict={}) -> None:
    """Runs the velocity analysis outputing a loom file

    BAMFILE or [BAMFILES] one or several bam files with position-sorted

    GTFFILE annotation file

    NOTE: it is keyword only argument function
    """
    
    ########################
    #    Resolve Inputs    #
    ########################

    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s',
                        level=[logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG][verbose])

    if isinstance(bamfile, tuple) and len(bamfile) > 1 and bamfile[-1][-4:] in [".bam", ".sam"]:
        multi = True
    elif isinstance(bamfile, tuple) and len(bamfile) == 1:
        multi = False
    else:
        raise IOError(f"Something went wrong in the argument parsing. You passed as bamfile: {bamfile}")

    if onefilepercell and multi:
        if bcfile is not None:
            raise ValueError("Inputs incompatibility. --bcfile/-b option was used together with --onefilepercell/-c option.")
        logging.warning("Each bam file will be interpreted as a DIFFERENT cell")
    elif not onefilepercell and multi:
        logging.warning("Several input files but --onefilepercell is False. Each bam file will be interpreted as containing a SET of cells!!!")

    if sampleid is None:
        assert metadatatable is None, "--metadatatable was specified but cannot fetch sample metadata without valid sampleid"
        if multi:
            logging.warning(f"When using mutliple files you may want to use --sampleid option to specify the name of the output file")
        if multi and not onefilepercell:
            full_name = "_".join([os.path.basename(bamfile[i]).split(".")[0] for i in range(len(bamfile))])
            if len(full_name) > 50:
                sampleid = f'multi_input_{os.path.basename(bamfile[0]).split(".")[0]}_{id_generator(5)}'
            else:
                sampleid = f'multi_input_{full_name}_and_others_{id_generator(5)}'
        elif multi and onefilepercell:
            sampleid = f'onefilepercell_{os.path.basename(bamfile[0]).split(".")[0]}_and_others_{id_generator(5)}'
        else:
            sampleid = f'{os.path.basename(bamfile[0]).split(".")[0]}_{id_generator(5)}'
        logging.info(f"No SAMPLEID specified, the sample will be called {sampleid} (last 5 digits are a random-id to avoid overwriting some other file by mistake)")

    # Create an output folder inside the cell ranger output folder
    if outputfolder is None:
        outputfolder = os.path.join(os.path.split(bamfile[0])[0], "velocyto")
        logging.info(f"No OUTPUTFOLDER specified, find output files inside {outputfolder}")
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    logic_class = getattr(vcy, logic)
    if not issubclass(logic_class, vcy.Logic):
        raise ValueError(f"{logic} is not a valid logic. Choose one among {', '.join([k for k, v in vcy.logic.__dict__.items() if issubclass(v, vcy.Logic)])}")
    else:
        logging.debug(f"Using logic: {logic}")
        logic_obj = logic_class()

    if bcfile is None:
        logging.debug("Cell barcodes will be determined while reading the .bam file")
        valid_bcset = None
    else:
        # Get valid cell barcodes
        valid_bcs_list = open(bcfile).read().rstrip().split()
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])  # with sample id and with -1
        if len(set(bc.split('-')[0] for bc in valid_bcs_list)) == 1:
            gem_grp = f"-{valid_bcs_list[0].split('-')[-1]}"
        else:
            gem_grp = "x"
        valid_bcset = set(bc.split('-')[0] for bc in valid_bcs_list)  # without -1
        logging.info(f"Read {len(valid_bcs_list)} cell barcodes from {bcfile}")
        logging.debug(f"Example of barcode: {valid_bcs_list[0].split('-')[0]} and cell_id: {valid_cellid_list[0]}")
        
    # Get metadata from sample sheet
    if metadatatable:
        try:
            sample_metadata = vcy.MetadataCollection(metadatatable)
            sample = sample_metadata.where("SampleID", sampleid)
            if len(sample) == 0:
                logging.error(f"Sample ID {sampleid} not found in sample sheet")
                # schema = []  # type: List
                sample = {}
            elif len(sample) > 1:
                logging.error(f"Sample ID {sampleid} has multiple lines in sample sheet")
                sys.exit(1)
            else:
                # schema = sample[0].types
                sample = sample[0].dict
            logging.debug(f"Collecting column attributes from {metadatatable}")
        except (NameError, TypeError) as e:
            logging.warn("SAMPLEFILE was not specified. add -s SAMPLEFILE to add metadata.")
            sample = {}
    else:
        sample = {}

    ########################
    #     Start Analysis   #
    ########################

    # Initialize Exon-Intron Counter with the logic and valid barcodes (need to do it now to peek)
    if without_umi:
        if umi_extension != "no":
            logging.warning("--umi-extension was specified but uncompatible with --without-umi, it will be ignored!")
        umi_extension = "without_umi"
    exincounter = vcy.ExInCounter(sampleid=sampleid, logic=logic_class, valid_bcset=valid_bcset, umi_extension=umi_extension, onefilepercell=onefilepercell, dump_option=dump, outputfolder=outputfolder)

    # Heuristic to chose the memory/cpu effort
    try:
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
    except subprocess.CalledProcessError:
        logging.warning("Your system does not support calling `grep MemAvailable /proc/meminfo` so the memory effort for the samtools command could not be chosen appropriatelly. 32Gb will be assumed")
        mb_available = 32000  # 64Gb
    threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
    mb_to_use = int(min(samtools_memory, mb_available / (len(bamfile) * threads_to_use)))
    compression = vcy.BAM_COMPRESSION

    # I need to peek into the bam file to know wich cell barcode flag should be used
    if onefilepercell:
        logging.debug("The multi input option ")
        tagname = "NOTAG"
    else:
        exincounter.peek(bamfile[0])
        tagname = exincounter.cellbarcode_str
    
    if multi and onefilepercell:
        bamfile_cellsorted = list(bamfile)
    else:
        bamfile_cellsorted = [f"{os.path.join(os.path.dirname(bmf), 'cellsorted_' + os.path.basename(bmf))}" for bmf in bamfile]

    sorting_process: Dict[int, Any] = {}
    for ni, bmf_cellsorted in enumerate(bamfile_cellsorted):
        # Start a subprocess that sorts the bam file
        command = f"samtools sort -l {compression} -m {mb_to_use}M -t {tagname} -O BAM -@ {threads_to_use} -o {bmf_cellsorted} {bamfile[ni]}"
        if os.path.exists(bmf_cellsorted):
            logging.warning(f"The file {bmf_cellsorted} already exists. The sorting step will be skipped and the existing file will be used.")
            check_end_process = False
        else:
            sorting_process[ni] = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
            logging.info(f"Starting the sorting process of {bamfile[ni]} the output will be at: {bmf_cellsorted}")
            logging.info(f"Command being run is: {command}")
            logging.info(f"While the bam sorting happens do other things...")
            check_end_process = True

    # Load annotations
    logging.info(f"Load the annotation from {gtffile}")
    annotations_by_chrm_strand = exincounter.read_transcriptmodels(gtffile)
    chrs = list(v for k, v in annotations_by_chrm_strand.items())
    tms = list(itertools.chain.from_iterable((v.values() for v in chrs)))
    ivls = list(itertools.chain.from_iterable(tms))
    logging.debug(f"Generated {len(ivls)} features corresponding to {len(tms)} transcript models from {gtffile}")
    del chrs, tms, ivls

    # Load annotations
    if repmask is not None:
        logging.info(f"Load the repeat masking annotation from {repmask}")
        mask_ivls_by_chromstrand = exincounter.read_repeats(repmask)

    # Go through the bam files a first time to markup introns
    logging.info(f"Scan {' '.join(bamfile)} to validate intron intervals")
    if test:  # NOTE: Remove this after finishing testing, the only purpuso was to save 15min in the debugging process
        logging.warning("This place is for developer only!")
        import pickle
        if os.path.exists("exincounter_dump.pickle"):
            logging.debug("exincounter_dump.pickle is being loaded")
            exincounter = pickle.load(open("exincounter_dump.pickle", "rb"))
        else:
            logging.debug("exincounter_dump.pickle was not found")
            logging.debug("Dumping exincounter_dump.pickle BEFORE markup")
            pickle.dump(exincounter, open("exincounter_dump.pickle", "wb"))
            exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)
    else:
        exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)

    # Wait for child process to terminate
    if check_end_process:
        logging.info(f"Now just waiting that the bam sorting process terminates")
        for k in sorting_process.keys():
            returncode = sorting_process[k].wait()
            if returncode == 0:
                logging.info(f"bam file #{k} has been sorted")
            else:
                raise MemoryError(f"bam file #{k} could not be sorted by cells.\n\
                This is probably related to an old version of samtools, please install samtools >= 1.6.\
                In alternative this could be a memory error, try to set the --samtools_memory option to a value compatible with your system. \
                Otherwise sort manually by samtools ``sort -l [compression] -m [mb_to_use]M -t [tagname] -O BAM -@ [threads_to_use] -o cellsorted_[bamfile] [bamfile]``")

    # Do the actual counting
    logging.debug("Start molecule counting!")
    results = exincounter.count(bamfile_cellsorted, multimap=multimap)  # NOTE: we would avoid some millions of if statements evalution if we write two function count and count_with output
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
        logging.debug(f"Example of barcode: {valid_bcs_list[0]} and cell_id: {valid_cellid_list[0]}")
     
    ca = {"CellID": np.array([f"{sampleid}:{v_bc}{gem_grp}" for v_bc in cell_bcs_order])}
    ca.update(additional_ca)

    for key, value in sample.items():
        ca[key] = np.full(len(cell_bcs_order), value)

    # Save to loom file
    outfile = os.path.join(outputfolder, f"{sampleid}.loom")
    logging.debug(f"Generating output file {outfile}")
    
    # row attributes
    atr_table = (("Gene", "genename", str),
                 ("Accession", "geneid", str),
                 ("Chromosome", "chrom", str),
                 ("Strand", "strand", str),
                 ("Start", "start", int),
                 ("End", "end", int))

    logging.debug("Collecting row attributes")
    ra = {}
    for name_col_attr, name_obj_attr, dtyp in atr_table:
        tmp_array = np.zeros((len(exincounter.genes),), dtype=object)  # type: np.ndarray
        for gene_id, gene_info in exincounter.genes.items():
            tmp_array[exincounter.geneid2ix[gene_id]] = getattr(gene_info, name_obj_attr)
        ra[name_col_attr] = tmp_array.astype(dtyp)
    
    logging.debug("Generating data table")
    layers: Dict[str, np.ndarray] = {}

    for layer_name in logic_obj.layers:
        layers[layer_name] = np.concatenate(dict_list_arrays[layer_name], axis=1)
        del dict_list_arrays[layer_name]
    
    for layer_name in logic_obj.layers:
        total: np.ndarray  # This is just a type annotation to avoid mypy compaints
        try:
            total += layers[layer_name]
        except NameError:
            total = np.array(layers[layer_name])

    logging.debug("Writing loom file")
    try:
        ds = loompy.create(filename=outfile, matrix=total, row_attrs=ra, col_attrs=ca, dtype="float32")
        for layer_name in logic_obj.layers:
            ds.set_layer(name=layer_name, matrix=layers[layer_name], dtype=vcy.LOOM_NUMERIC_DTYPE)
        ds.attrs["velocyto.__version__"] = vcy.__version__
        ds.attrs["velocyto.logic"] = logic
        ds.close()
    except TypeError:
        # If user is using loompy2
        # NOTE maybe this is not super efficient if the type and order are already correct
        tmp_layers = {"": total.astype("float32", order="C", copy=False)}
        tmp_layers.update({layer_name: layers[layer_name].astype(vcy.LOOM_NUMERIC_DTYPE, order="C", copy=False) for layer_name in logic_obj.layers})
        loompy.create(filename=outfile, layers=tmp_layers, row_attrs=ra, col_attrs=ca, file_attrs={"velocyto.__version__": vcy.__version__, "velocyto.logic": logic})
    logging.debug("Terminated Succesfully!")
