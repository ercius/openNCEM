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


Requirements
------------

``ncempy`` is designed and written for python3.5.

It relies on the following packages:
* numpy
* scipy
* matplotlib (for plotting)
* h5py (for EMD files)
* h5py_cache (for EMD Velox files)

Installation
------------

We recommend installing the ``ncempy`` package from PyPi:

``pip install ncempy``

License
-------

``ncempy`` is dual licensed under GPLv3 and MIT. The io module is the only part
released under MIT to improve interoperability with other packages.
