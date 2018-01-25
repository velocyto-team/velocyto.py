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

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


def id_generator(size: int=6, chars: str=string.ascii_uppercase + string.digits) -> str:
    return ''.join(random.choice(chars) for _ in range(size))


def _run(*, bamfile: str, gtffile: str,
         bcfile: str, outputfolder: str,
         sampleid: str, metadatatable: str,
         repmask: str, logic: str, molrep: bool,
         multimap: bool, test: bool, samtools_threads: int, samtools_memory: int,
         additional_ca: dict={}) -> None:
    """Runs the velocity analysis outputing a loom file

    BAMFILE bam file with sorted reads

    GTFFILE annotation file

    NOTE: it is keyword only argument function
    """
    
    ########################
    #    Resolve Inputs    #
    ########################

    if sampleid is None:
        assert metadatatable is None, "Cannot fetch sample metadata without valid sampleid"
        sampleid = f'{os.path.basename(bamfile).split(".")[0]}_{id_generator(5)}'
        logging.debug(f"No SAMPLEID specified, the sample will be called {sampleid}")

    # Create an output folder inside the cell ranger output folder
    if outputfolder is None:
        outputfolder = os.path.join(os.path.split(bamfile)[0], "velocyto")
        logging.debug(f"No OUTPUTFOLDER specified, find output files inside {outputfolder}")
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    logic_obj = getattr(vcy, logic)
    if not issubclass(logic_obj, vcy.Logic):
        raise ValueError(f"{logic} is not a valid logic. Chose one among {', '.join([k for k, v in vcy.logic.__dict__.items() if issubclass(v, vcy.Logic)])}")
    else:
        logging.debug(f"Using logic: {logic}")

    if bcfile is None:
        logging.debug("Cell barcodes will be determined while reading the .bam file")
        valid_bcset = None
    else:
        # Get valid cell barcodes
        valid_bcs_list = [l.strip() for l in open(bcfile).readlines()]
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])  # with sample id and with -1
        if len(set(bc.split('-')[0] for bc in valid_bcs_list)) == 1:
            gem_grp = f"-{valid_bcs_list[0].split('-')[-1]}"
        else:
            gem_grp = "x"
        valid_bcset = set(bc.split('-')[0] for bc in valid_bcs_list)  # without -1
        logging.debug(f"Read {len(valid_bcs_list)} cell barcodes from {bcfile}")
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
    exincounter = vcy.ExInCounter(logic_obj, valid_bcset)

    # Heuristic to chose the memory/cpu effort
    mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
    threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
    mb_to_use = int(min(samtools_memory, mb_available / threads_to_use))
    compression = vcy.BAM_COMPRESSION

    # I need to peek into the bam file to know wich cell barcode flag should be used
    exincounter.peek(bamfile)
    tagname = exincounter.cellbarcode_str
    bamfile_cellsorted = f"{os.path.join(os.path.dirname(bamfile), 'cellsorted_' + os.path.basename(bamfile))}"

    # Start a subprocess that sorts the bam file
    command = f"samtools sort -l {compression} -m {mb_to_use}M -t {tagname} -O BAM -@ {threads_to_use} -o {bamfile_cellsorted} {bamfile}"
    if os.path.exists(bamfile_cellsorted):
        logging.warning(f"The file {bamfile_cellsorted} already exists. The sorting step will be skipped and the existing file will be used.")
        check_end_process = False
    else:
        sorting_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        logging.info(f"Starting the sorting process of {bamfile} the output will be at: {bamfile_cellsorted}")
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

    # Go through the sam files a first time to markup introns
    logging.info(f"Scan {bamfile} to validate intron intervals")
    if test:  # NOTE: Remove this after finishing testing, the only purpuso was to save 15min in the debugging process
        import pickle
        if os.path.exists("exincounter_dump.pickle"):
            logging.debug("exincounter_dump.pickle is being loaded")
            exincounter = pickle.load(open("exincounter_dump.pickle", "rb"))
        else:
            logging.debug("exincounter_dump.pickle was not found")
            logging.debug("Dumping exincounter_dump.pickle BEFORE markup")
            pickle.dump(exincounter, open("exincounter_dump.pickle", "wb"))
            exincounter.mark_up_introns(samfile=bamfile, multimap=multimap)
            
    else:
        exincounter.mark_up_introns(samfile=bamfile, multimap=multimap)

    # Wait for child process to terminate
    if check_end_process:
        logging.info(f"Now just waiting that the bam sorting process terminates")
        sorting_process.wait()
        logging.info(f"bam file has been sorted")

    # Do the actual counting
    logging.debug("Start molecule counting!")
    results = exincounter.count(bamfile_cellsorted, multimap=multimap, molecules_report=molrep)  # NOTE: we would avoid some millions of if statements evalution if we write two function count and count_with output
    list_spliced_arrays, list_unspliced_arrays, list_ambiguous_arrays, cell_bcs_order = results

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
        ca[key] = np.array([value] * len(cell_bcs_order))

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
    spliced = np.concatenate(list_spliced_arrays, axis=1)
    del list_spliced_arrays
    unspliced = np.concatenate(list_unspliced_arrays, axis=1)
    del list_unspliced_arrays
    ambiguous = np.concatenate(list_ambiguous_arrays, axis=1)
    del list_ambiguous_arrays
    
    total = spliced + unspliced + ambiguous
    logging.debug("Writing loom file")
    try:
        ds = loompy.create(filename=outfile, matrix=total, row_attrs=ra, col_attrs=ca, dtype="float32")
        ds.set_layer(name="spliced", matrix=spliced, dtype=vcy.LOOM_NUMERIC_DTYPE)
        ds.set_layer(name="unspliced", matrix=unspliced, dtype=vcy.LOOM_NUMERIC_DTYPE)
        ds.set_layer(name="ambiguous", matrix=ambiguous, dtype=vcy.LOOM_NUMERIC_DTYPE)
        ds.attrs["velocyto.__version__"] = vcy.__version__
        ds.close()
    except TypeError:
        # If user is using loompy2
        ds = loompy.create(filename=outfile, layers={"": total.astype("float32", order="C", copy=False),
                                                     "spliced": spliced.astype(vcy.LOOM_NUMERIC_DTYPE, order="C", copy=False),
                                                     "unspliced": unspliced.astype(vcy.LOOM_NUMERIC_DTYPE, order="C", copy=False),
                                                     "ambiguous": ambiguous.astype(vcy.LOOM_NUMERIC_DTYPE, order="C", copy=False)},
                           row_attrs=ra, col_attrs=ca, file_attrs={"velocyto.__version__":vcy.__version__})
    logging.debug("Terminated Succesfully!")
