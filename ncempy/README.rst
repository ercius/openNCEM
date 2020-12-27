------
ncempy
------

openNCEM's python package.

Structure
---------

The package is given the following structure:

* algo
    Algorithms used for image processing and other computing. These act as the machinery of the provided tools. Heavy reusing is encouraged by keeping them general to the processing of datasets.

* eval
    Evaluation routines build from the single algorithms in `algo`. These address specific tasks like evaluating the results from a particular method or experimental setup.

* edstomo
    (Optional) Routines used for reconstructing STEM/EDS tomographic datasets.

* io
    Module to do file IO for various file formats. While the EMD file format is used internally, other file formats commonly used in electron microscopy are read in using importers.

* data
    Example datasets which can be used to test and demonstrate the code base.

* docs
    Documentation of the `ncempy` package.

* test
    Tests for all modules, functions, and classes.


Requirements
------------

``ncempy`` is designed and written for python3.6 or later.

It relies on the following packages:
* numpy
* scipy
* matplotlib (for plotting)
* h5py (for EMD files)

edstomo has additional optional packages:
* glob2
* genfire
* hyperspy
* scipy
* scikit-image
* matplotlib
* ipyvolume

Installation
------------

We recommend installing the ``ncempy`` package from PyPi:

``pip install ncempy``

If you wish to install the optional EDSTomo module then run:

``pip install 'ncempy[edstomo]'``

License
-------

``ncempy`` is dual licensed under GPLv3 and MIT. The io module is the only part
released under MIT to improve interoperability with other packages.
