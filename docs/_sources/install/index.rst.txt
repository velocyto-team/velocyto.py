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

To make sure you have all the dependencies correctly installed, including libraries trickier to install do:

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click loompy

.. _pypi:

Install using PyPI
------------------

Just run:

::

    pip install velocyto

`pysam` will be installed by ``pip`` as dependencies

.. tip::
    ``velocyto 0.9.x`` is an alpha release and it is updated often. If you installed with pip make sure you run ``pip install -U --no-deps velocyto`` often.

.. _command:

Install from source
~~~~~~~~~~~~~~~~~~~

After making sure all the requirements are satisfied. You can just run:

::

    git clone https://github.com/velocyto-team/velocyto.py.git
    cd velocyto.py
    pip install -e .  # installs in development mode, making symlinks to the current directory


This command will install the python library ``velocyto`` and the homonymous command line tool.

At this point you can test the installation was successful by running ``velocyto --help``.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 

.. tip::
    ``velocyto 0.9.x`` is an alpha release, we recommend to pull the latest bufixes and feature improvements often. Adding the ``-e`` flags installs the software in `development` mode, you can simple do ``git pull`` to update the library.


.. _conda:

Install using conda
-------------------

.. note::
   This installation method is not currently available. The plan is make it available upon 1.0 release.
