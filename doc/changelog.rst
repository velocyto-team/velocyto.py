.. _changelog:

=========
Changelog
=========
* :release:`0.13.1 <2018-01-22>`
* :feature:`-` Sample metadata file can be specified with different csv formats (the format will be determined automatically)
* :release:`0.13.0 <2018-01-19>`
* :bug:`-` Sometimes velocyto missed to detect and warn the user that the `.gtf` genome annotation file was not sorted, this could have caused undetected errors in the analysis. If you run velocyto without sorting the .gtf, we suggest rerunning.
* :feature:`-` CLI does not require presorting the gtf files. To reduce possibility of incorrect usage, now .gtf file sorting sorting is performed in memory (and not saved).
* :feature:`-` Improve documentation: remove information about sorting .gtf files. This procedure is not needed anymore.
* :release:`0.12.4 <2018-01-18>`
* :bug:`40` Error in hdf5 serialization when using cluster label as object array is now fixed
* :release:`0.12.3 <2018-01-17>`
* :bug:`38` Error in running run10x is now fixed
* :release:`0.12.2 <2018-01-12>`
* :bug:`37` Initial cell size array gets updated properly when filtering cells
* :release:`0.12.1 <2018-01-04>`
* :bug:`35` Attempted to fix a reported bug running velocyto CLI
* :release:`0.12.0 <2017-12-17>`
* :feature:`-` Add possibility to export pickle containing information of every molecule
* :feature:`-` Remove the subcommand ``multi10x``
* :bug:`- major` Incorrect 0-based indexing for splicing junction corrected (was not causing problems because buffered by MIN_FLANK) 
* :bug:`- major` Many small bug fixes
* :bug:`31 major` Memory usage bug should be solved.
* :feature:`-` Large parts of the documentation rewritten to match the changes in API
* :feature:`-` Make the CLI simpler removing the extract interval step. 
  Now the source .gtf files can be provided directly, they should be provided sorted using :code:`sort -k1,1 -k7,7 -k4,4n -o [OUTFILE] [INFILE]`
* :feature:`-` Changelog added to the doc
* :support:`-` Update the documentation for the new  :ref:`CLI <cli>`
* :feature:`-` Support different Logic levels
* :feature:`-` Pipeline now consider all the possible transcript models that could be supported by a set of reads individually and then decides on the spliced/unspliced/ambigous count.
* :release:`0.11.0 <2017-12-01>`
* :bug:`- major` fix a bug with ambiguous molecules counting and version bump
* :release:`0.10.3 <2017-11-23>`
* :bug:`- major` The debug and sampleid option had the same short flag `-d`
* :release:`0.10.2 <2017-11-18>`
* :release:`0.10.1 <2017-11-18>`
* :feature:`-` further ~5x speedup of cython functions making them 100% C and using malloc instead of memory views
* :release:`0.10.0 <2017-11-18>`
* :feature:`-` Add support for DropSeq pipelines where the barcode flags in the bam file are `XC` and `XM` instead of `CB` and `UB`
* :bug:`- major` Using sphinx 1.7 sorts the autodoc API correctly
* :release:`0.9.13 <2017-11-04>`
* :release:`0.9.12 <2017-11-04>`
* :release:`0.9.11 <2017-11-03>`
* :feature:`-` Improve the docs
* :release:`0.9.10 <2017-11-02>`
* :release:`0.9.9 <2017-10-31>`
* :release:`0.9.8 <2017-10-26>`
* :release:`0.9.7 <2017-10-25>`
* :release:`0.9.6 <2017-10-25>`