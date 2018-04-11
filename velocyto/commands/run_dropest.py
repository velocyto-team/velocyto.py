import sys
import os
import glob
import re
import gzip
import click
import loompy
import numpy as np
import random
import string
import csv
from collections import defaultdict
import logging
from typing import *
import velocyto as vcy
from ._run import _run

# logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


def id_generator(size: int=6, chars: str=string.ascii_uppercase + string.digits) -> str:
    return ''.join(random.choice(chars) for _ in range(size))


@click.command(short_help="Runs the velocity analysis on DropEst preprocessed data")
@click.argument("bamfile", nargs=1, required=True,
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("gtffile", required=True,
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.option("--bcfile", "-b",
              help="""Valid barcodes file, to filter the bam. If --bcfile is not specified the file will be searched in the default position outputed by ``velocyto tools dropest_bc_correct``. Otherwise an error will be thrown""",
              default=None,
              show_default=True,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--logic", "-l",
              help="The logic to use for the filtering (default: Default)",
              default="Default")
@click.option("--outputfolder", "-o",
              help="Output folder, if it does not exist it will be created.",
              default=None,
              type=click.Path(exists=False))
@click.option("--sampleid", "-e",
              help="The sample name that will be used as a the filename of the output.",
              default=None,
              type=click.Path(exists=False))
@click.option("--repmask", "-m",
              help=".gtf file containing intervals to mask (Optional)",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--samtools-threads", "-@",
              help="The number of threads to use to sort the bam by cellID file using samtools",
              default=16)
@click.option("--samtools-memory",
              help="The number of MB used for every thread by samtools to sort the bam file",
              default=2048)
@click.option("--dump", "-d",
              help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed (default: 0)",
              default="0")
@click.option('--verbose', '-v',
              help="Set the vebosity level: -v (only warinings) -vv (warinings and info) -vvv (warinings, info and debug)",
              count=True, default=1)
def run_dropest(bamfile: str, gtffile: str, bcfile: str, logic: str, outputfolder: str, sampleid: str,
                repmask: str, samtools_threads: int, samtools_memory: int, dump: str, verbose: int, additional_ca: dict={}) -> None:
    """Runs the velocity analysis on DropEst preprocessed data

    BAMFILE  bam files to be analyzed

    GTFFILE genome annotation file
    """
    if bcfile is None:
        parentpath, bamfilename = os.path.split(bamfile)
        bcfile = os.path.join(parentpath, f"barcodes_{bamfilename.split('_')[0]}.tsv")
        logging.info(f"Attempting to find automatically the valid barcode list file {bcfile}")
        if os.path.exists(bcfile):
            logging.info(f"{bcfile} found ")
            pass
        else:
            logging.info(f"{bcfile} not found!")
            logging.error("In ``run_dropest`` specifying the ``--bcfile/-b`` is required. Use ``run`` if you want to adventure in a more custom usage.")
            logging.info("Exit without doing nothing")
            return

    if "correct" not in bamfile:
        logging.warning("The file you are using does not start with the prefix ``correct_`` so it might not be the output of ``velocyto tools dropest_bc_correct``.")
        logging.info("The program will run despite the warning but be aware of the possible consequences of not correcting the barcodes")
    return _run(bamfile=bamfile, gtffile=gtffile, bcfile=bcfile, outputfolder=outputfolder,
                sampleid=sampleid, metadatatable=None, repmask=repmask, onefilepercell=False,
                logic=logic, without_umi=False, umi_extension="chr",
                multimap=False, test=False, samtools_threads=samtools_threads,
                samtools_memory=samtools_memory, dump=dump, verbose=verbose, additional_ca=additional_ca)
