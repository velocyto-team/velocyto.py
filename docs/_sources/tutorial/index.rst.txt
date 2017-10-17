.. _tutorial:

Tutorial
========

Velocyto consists of two main components:

    - A command line interface (CLI), that is used to run the pipeline that generates spliced/unspliced expression matrices.
    
    - A library including functions to estimate RNA velocity from the above mentioned data matrices.

Getting Started
---------------
First make sure that velocyto is correctly installed following :ref:`this guide <install>`.


Running the CLI
---------------
This guide explains how to run the counting pipeline.
This starts from .bam/.sam files and yields a `.loom file <https://github.com/linnarsson-lab/loompy>`_ with counts divided in the cathegories spliced/unspliced/ambiguous. 

.. toctree::
    cli

Estimating RNA velocity
-----------------------
This guide covers the analysis and assumes that you have produced a `.loom file <https://github.com/linnarsson-lab/loompy>`_ following the guide above.

.. toctree::
    analysis

