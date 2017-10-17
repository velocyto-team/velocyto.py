.. _install:

Installation Guide
==================

.. _fromsource:

Install from source
-------------------

.. _require:

Requirements
~~~~~~~~~~~~

To run velocyto you will need python >=3.6.0 (we have no plans to support python<=3.5).
We recommend to use `anaconda <https://www.continuum.io/downloads>`_ and the ``conda`` command to install dependencies limiting the use of ``pip`` to libraries not available on conda's main channel. (``pip`` can be used too but its dependency-managing migh be less robust). 

To make sure you have all the dependencies first do:

::

    conda install numpy scipy cython scikit-learn scikit-image h5py click
    pip install pysam loompy

.. _command:

Install command
~~~~~~~~~~~~~~~

After all the requirements are satisfied. You can just run:

::

    git clone https://github.com/linnarsson-lab/velocyto.git
    cd velocyto/python
    pip install -e .


This command will install the python library ``velocyto`` and the homonymous command line tool.
The ``-e`` flags installs the software in `development` mode, so to update you can simple do ``git pull``.

At this point you can test the installation was successful by running ``velocyto --help``.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 

.. _pypi:

Install using PyPI wheels
-------------------------

.. note::
   This installation modality is not currently available. The plan is make it available upon 1.0 release.

.. _conda:

Install using conda
-------------------

.. note::
   This installation modality is not currently available. The plan is make it available upon 1.0 release.