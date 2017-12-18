import sys
import os
import subprocess
import logging
import click
import time
from typing import *


@click.command(short_help="Sorts one or more bam files in cell-barcode order.")
@click.argument("bamfiles", nargs=-1, type=click.Path(exists=True))
@click.option("--serial", "-s", default=False, is_flag=True, help="Serial implementation, the arguments regarding multiprocessign will be ignored.")
@click.option("--processes", "-p", type=int, default=3, help="Number of sorting to perform at the same time. (default: 3)")
@click.option("--threads", "-@", type=int, default=16, help="Number of threads to use for each sorting process. (default: 16)")
@click.option("--mb-per-thread", "-m", default=4096, help="Max MBs to use for each threads. (default: 4096)")
@click.option("--tagname", "-t", type=str, default="CB", help="Cell barcode tag used in the bam file. (default: CB)")
@click.option("--compression", "-l", type=int, default=7, help="compression of the bam file. (default: 7)")
def cellsort(bamfiles: Tuple, serial: bool, processes: int, threads: int, mb_per_thread: int, tagname: str, compression: str) -> None:
    """Sorts one or more bam files in cell-barcode order in a very fast and parallelized way.

    BAMFILES    one or a list of bamfiles

    Notes:

    It uses almost al the resources of the node to run several multithreaded processes.
    It makes a copy and never overwrites the original bam file.
    """
    file_list: List[str] = list(bamfiles)

    if serial:
        if type(bamfiles) is not tuple:
            raise NotImplementedError("multiple argument list expected")
        for bamfile in bamfiles:
            bamfile_cellsorted = f"{os.path.join(os.path.dirname(bamfile), 'cellsorted_' + os.path.basename(bamfile))}"
            mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
            mb_per_process = int(mb_available * 0.70)
            mb_per_thread = min(mb_per_thread, int(mb_per_process * 0.90 / threads))
            command = f"samtools sort -l {compression} -m {mb_per_thread}M -t {tagname} -O BAM -@ {threads} -o {bamfile_cellsorted} {bamfile}"
            if os.path.exists(bamfile_cellsorted):
                logging.warning(f"The file {bamfile_cellsorted} already exist. File will not be overwritten! Closing.")
            else:
                logging.info(f"Starting the sorting process of {bamfile} the output will be at: {bamfile_cellsorted}")
                logging.info(f"Command being run is: {command}")
                p = subprocess.check_call(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
                logging.info(f"Done")
    else:
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
        mb_per_process = int(mb_available * 0.80 / processes)
        mb_per_thread = min(mb_per_thread, int(mb_per_process * 0.90 / threads))

        logging.info(f"Launching {processes} processes with {threads} threads and {mb_per_thread}MB each.")
        sorting_processes: List = []
        while True:
            if len(sorting_processes) < processes and len(file_list) > 0:
                bamfile = file_list.pop(0)
                bamfile_cellsorted = f"{os.path.join(os.path.dirname(bamfile), 'cellsorted_' + os.path.basename(bamfile))}"
                command = f"samtools sort -l {compression} -m {mb_per_thread}M -t {tagname} -O BAM -@ {threads} -o {bamfile_cellsorted} {bamfile}"
                if os.path.exists(bamfile_cellsorted):
                    logging.warning(f"The file {bamfile_cellsorted} already exist. File will not be overwritten! Closing.")
                else:
                    sorting_processes.append(subprocess.Popen(command.split()))
                    logging.info(f"Starting the sorting process of {bamfile} the output will be at: {bamfile_cellsorted}")
                    logging.info(f"Command being run is: {command}")
                    logging.info(f"A total of {len(file_list)} files are waiting for processing")
            else:
                for i in range(len(sorting_processes)):
                    if sorting_processes[i].poll() is not None:  # not check_pid(sorting_processes[i].pid):
                        logging.info(f"The process PID:{sorting_processes[i].pid} is done.")
                        if check_pid(sorting_processes[i].pid):
                            logging.warning(f"PID:{sorting_processes[i].pid} still exists!")
                            break
                        # sorting_processes[i].kill()
                        del sorting_processes[i]
                        logging.info(f"A total of {len(file_list)} files are waiting for processing")
                        break

            if len(file_list) == 0 and len(sorting_processes) == 0:
                logging.info("Sorting finished")
                break
            else:
                time.sleep(3)


def check_pid(pid: int) -> bool:
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True
