.. _install:

Installation Guide
==================

Instructions for installing the python library ``velocyto`` and the command line tool of the same name.

.. _require:

Requirements
------------

To run velocyto you will need python >=3.6.0 (we have no plans to support python<=3.5).
We recommend to use `anaconda <https://www.continuum.io/downloads>`_ and the ``conda`` command to install dependencies limiting the use of ``pip`` to libraries not available on conda's main channel. (``pip`` can be used too but its dependency-managing might be less robust). 

To make sure you have all the dependencies correctly installed, including libraries trickier to install, run:

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click loompy


.. _pypi:

Install using PyPI
------------------

Just run:

::

    pip install velocyto

`pysam` will be installed by ``pip`` as a dependency

.. tip::
    ``velocyto 0.10.x`` is an alpha release and it is updated often. If you installed with pip make sure you run ``pip install -U --no-deps velocyto`` often.

You can test whether the installation was successful by running ``velocyto --help``.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 


.. _fromsource:

Install from source
-------------------

Make sure all the requirements are satisfied, and that `git` is installed on your system. After doing so, run the following commands to complete the installation:

::

    git clone https://github.com/velocyto-team/velocyto.py.git
    cd velocyto.py
    pip install -e .  # installs in development mode, making symlinks to the current directory

.. tip::
    ``velocyto 0.10.x`` is an alpha release, we recommend pulling in the latest bufixes and feature improvements often. Adding the ``-e`` flag installs the software in `development` mode, after which you can simply use ``git pull`` to update the library.

You can test whether the installation was successful by running ``velocyto --help``.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 


.. _conda:

Install using conda
-------------------

.. note::
   This installation method is not currently available. The plan is make it available upon the 1.0 release.
