ncempy package
==============

openNCEM's python package.

Structure
---------

The ncempy package bundles a set of subpackages:

* :doc:`ncempy.algo`
    Algorithms used for image processing and other computing. These act as the machinery of the provided tools. Heavy reusing is encouraged by keeping them general to the processing of datasets.

* :doc:`ncempy.eval`
    Evaluation routines build from the single algorithms in `algo`. These address specific tasks like evaluating the results from a particular method or experimental setup.

* :doc:`ncempy.io`
    Module to do file IO for various file formats. While the EMD file format is used internally, other file formats commonly used in electron microscopy are read in using importers.    

Requirements
------------

``ncempy`` is designed and written for python3.5.

It relies on the following packages:

    * numpy
    * scipy
    * matplotlib (for plotting)
    * h5py (for EMD files)

Installation
------------

For now we support pip installing the ``ncempy`` package from the gitHub repository:

``pip install 'git+https://github.com/ercius/openNCEM.git@development#egg=ncempy'``

Indices and tables
------------------
        
* :ref:`genindex`
* :ref:`modindex`

Contents
--------

.. toctree::
    :maxdepth: 1
    
    ncempy.algo
    ncempy.eval
    ncempy.io

