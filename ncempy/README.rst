------
ncempy
------

openNCEM's python package.

Structure
---------

The package is given the following structure:

* docs
    Documentation of the `ncempy` package.

* algo
    Algorithms used for image processing and other computing. These act as the machinery of the provided tools. Heavy reusing is encouraged by keeping them general to the processing of datasets.

* eval
    Evaluation routines build from the single algorithms in `algo`. These address specific tasks like evaluating the results from a particular method or experimental setup.

* io
    Module to do file IO for various file formats. While the EMD file format is used internally, other file formats commonly used in electron microscopy are read in using importers.

* test
    Unittests for all modules, functions, lines of code.

The documentation of ``ncempy`` can be found at https://openncem.readthedocs.io/en/latest/.

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

