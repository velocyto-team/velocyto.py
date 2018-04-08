import sys
import os
import glob
import re
import click
import numpy as np
import random
import string
import logging
from typing import *
import velocyto as vcy
from ._run import _run

# logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


def id_generator(size: int=6, chars: str=string.ascii_uppercase + string.digits) -> str:
    return ''.join(random.choice(chars) for _ in range(size))


@click.command(short_help="Runs the velocity analysis outputing a loom file")
@click.argument("bamfile", nargs=-1, required=True,
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
@click.option("--bcfile", "-b",
              help="""Valid barcodes file, to filter the bam. If --bcfile is not specified all the cell barcodes will be incuded.
              Cell barcodes should be specified in the bcfile as the `CB` tag for each read""",
              default=None,
              show_default=True,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--outputfolder", "-o",
              help="Output folder, if it does not exist it will be created.",
              default=None,
              type=click.Path(exists=False))
@click.option("--sampleid", "-e",
              help="The sample name that will be used to retrieve informations from metadatatable",
              default=None,
              type=click.Path(exists=False))
@click.option("--metadatatable", "-s",
              help="Table containing metadata of the various samples (csv formatted, rows are samples and cols are entries)",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--mask", "-m",
              help=".gtf file containing intervals to mask",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--onefilepercell", "-c",
              help="""If this flag is used every bamfile passed is interpreted as an independent cell, otherwise multiple files are interpreted as batch of different cells to be analized together.
              Important: cells reads should not be distributed over multiple bamfiles is not supported!! (default: off)""",
              default=False,
              is_flag=True)
@click.option("--logic", "-l",
              help="The logic to use for the filtering (default: Default)",
              default="Default")
@click.option("--without-umi", "-U",
              help="If this flag is used the data is assumed UMI-less and reads are counted instead of molecules (default: off)",
              default=False,
              is_flag=True)
@click.option("--umi-extension", "-u",
              help="""In case UMI is too short to guarantee uniqueness (without information from the ampping) set this parameter to `chr`, `Gene` ro `[N]bp`
              If set to `chr` the mapping position (binned to 10Gb intervals) will be appended to `UB` (ideal for InDrops+dropEst). If set to `Gene` then the `GX` tag will be appended to the `UB` tag.
              If set to `[N]bp` the first N bases of the sequence will be used to extend `UB` (ideal for STRT). (Default: `no`)""",
              default="no")
@click.option("--multimap", "-M",
              help="""Consider not unique mappings (not reccomended)""",
              default=False,
              is_flag=True)
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
def run(bamfile: str, gtffile: str,
        bcfile: str, outputfolder: str,
        sampleid: str, metadatatable: str,
        mask: str, onefilepercell: bool, logic: str, without_umi: str, umi_extension: str, multimap: bool,
        samtools_threads: int, samtools_memory: int, dump: str, verbose: int,
        additional_ca: dict={}) -> None:
    """Runs the velocity analysis outputing a loom file

    BAMFILE bam file with sorted reads

    GTFFILE genome annotation file
    """
    return _run(bamfile=bamfile, gtffile=gtffile, bcfile=bcfile, outputfolder=outputfolder,
                sampleid=sampleid, metadatatable=metadatatable, repmask=mask, onefilepercell=onefilepercell,
                logic=logic, without_umi=without_umi, umi_extension=umi_extension, multimap=multimap, test=False, samtools_threads=samtools_threads,
                samtools_memory=samtools_memory, dump=dump, verbose=verbose, additional_ca=additional_ca)
