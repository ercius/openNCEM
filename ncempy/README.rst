======
ncempy
======

openNCEM's python package

=========
Structure
=========

The package is given the following structure:

* docs

 Documentation of the `ncempy` package.

* algo

 Algorithms used for image processing and other computing. These act as the machinery of the provided tools. Heavy reusing is encouraged by keeping them general to the processing of datasets.

* eval

 Evaluation routines build from the single algorithms in `algo`. These address specific tasks like evaluating the results from a particular method or experimental setup.

* fio

 Module to do file IO for various file formats. While the EMD file format is used internally, other file formats commonly used in electron microscopy are read in using importers.


* test

 Unittests for all modules, functions, lines of code.

See the readmes in the subfolders for details.

============
Requirements
============

``ncempy`` is designed for python3.x.

It relies on the following packages:

* numpy
* scipy
* matplotlib (for plotting)
* PySide (for GUI)
* pygtgraph (for plotting in GUI)
* h5py (for EMD files)

============
Installation
============

For now we support pip installing it from the gitHub repository:

``pip install 'git+https://github.com/ercius/openNCEM.git@development#egg=ncempy&subdirectory=ncempy'``


Problem with PySide not being available for python3.5, somehow it works on ubuntu if you use it from a system package. So I create a virtualenv with access to system site packages:

``mkvirtualenv ncempy -p python3 --system-site-packages``