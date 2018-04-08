.. _install:

Installation Guide
==================

Please follow this guide to install both the ``velocyto``  command line tool and the homonymous library for interactive analysis.

.. _require:

Requirements
------------

To run velocyto you will need python >=3.6.0 (we have no plans to support python<=3.5).
We recommend using `anaconda <https://www.continuum.io/downloads>`_ and the ``conda`` command to install dependencies (of course you can use ``pip`` but be aware its dependency-managing might be less robust).
Feel free to use ``pip`` if some libraries are not available on the conda channels you are using. 

Run the following command to make sure you have all the dependencies correctly installed:

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click

Note: ``pysam`` is also a requirement but currently it is preferable to install it through PyPI wheel (i.e. by ``pip install pysam``) as we have experienced problems with the version packaged by conda-forge.

.. _pypi:

Install using PyPI
------------------

Every new stable version of ``velocyto`` gets immediately released on PyPI, so running the following command will install on your system the cutting-edge version:

::

    pip install velocyto

`pysam` and `loompy` will be installed by ``pip`` as a dependencies.

.. tip::
    `velocyto` |version| is an alpha release and it is updated often. If you installed with pip make sure you run ``pip install -U --no-deps velocyto`` now and then.

You can test whether the installation was successful by running ``velocyto --help``.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 


.. _fromsource:

Install from source
-------------------

If you plan to explore and make changes to the source code, or you have requested some bug-fix that is temporarily available only on the github ``dev`` branch, then you need to install ``velocyto`` directly from source.


First of all, make sure all the dependencies are installed, and that `git` is installed on your system. 
Then, run the following commands to complete the installation:

::

    git clone https://github.com/velocyto-team/velocyto.py.git
    cd velocyto.py
    pip install -e .  # note the trailing dot

You can test whether the installation was successful by running ``velocyto --help``.

.. tip::
    `velocyto` |version| is an alpha release, we recommend pulling in the latest bufixes and feature improvements often. Note that adding the ``-e`` flag to the pip command installs the software in `development` mode, when a package is installed this way each change to the source immediatelly reflects to changes in the installed library. This means that you can simply use ``git pull`` to update your installation.

To get started with ``velocyto`` you can follow :ref:`our guide <tutorial>`. 


.. _conda:

Install using conda
-------------------

.. note::
   This installation method is not currently available. Our plan is make it available upon the 1.0 release.


Extra requirements
------------------

If you want to use velocyto in combination with DropEst you will have to install some extra requirements: ``R`` and ``rpy2``.
We suggest to this with ``conda``:

::

    conda install R rpy2


