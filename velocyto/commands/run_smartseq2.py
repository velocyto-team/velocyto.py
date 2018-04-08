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


@click.command(short_help="Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)")
@click.argument("bamfiles", nargs=-1, required=True,
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.option("--outputfolder", "-o",
              help="Output folder, if it does not exist it will be created.",
              default=None,
              type=click.Path(exists=False))
@click.option("--sampleid", "-e",
              help="The sample name that will be used as a the filename of the output.",
              default=None,
              type=click.Path(exists=False))
@click.option("--repmask", "-m",
              help=".gtf file containing intervals to mask",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--dump", "-d",
              help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed (default: 0)",
              default="0")
@click.option('--verbose', '-v',
              help="Set the vebosity level: -v (only warinings) -vv (warinings and info) -vvv (warinings, info and debug)",
              count=True, default=1)
def run_smartseq2(bamfiles: str, gtffile: str, outputfolder: str, sampleid: str,
                  repmask: str, dump: str, verbose: int, additional_ca: dict={}) -> None:
    """Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)

    [BAMFILES, ...] a sequence of bam files to be analyzed (e.g. use a wild-card expansion)

    GTFFILE genome annotation file
    """
    return _run(bamfile=bamfiles, gtffile=gtffile, bcfile=None, outputfolder=outputfolder,
                sampleid=sampleid, metadatatable=None, repmask=repmask, onefilepercell=True,
                logic="SmartSeq2", without_umi=True, umi_extension="no",
                multimap=False, test=False, samtools_threads=1,
                samtools_memory=1, dump=dump, verbose=verbose, additional_ca=additional_ca)
