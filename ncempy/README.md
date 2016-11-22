# ncempy

openNCEM's python package


## Structure

The package is given the following structure:

- docs


- algo

 Algorithms used for image processing and other computing. These act as the machinery of the provided tools. Heavy reusing is encouraged by keeping them general to the processing of datasets.

- eval

 Evaluation routines build from the single algorithms in `algo`. These address specific tasks like evaluating the results from a particular method or experimental setup.

- io

 Module to do file IO for various file formats. While the EMD file format is used internally, other file formats commonly used in electron microscopy are read in using importers.


- test

 Unittests for all modules, functions, lines of code.

See the readmes in the subfolders for details.


## Requirements

`ncempy` is designed for python3.x.

It relies on the following packages:
- `numpy`
- `scipy`
- `matplotlib` (for plotting)
- `PySide` (for GUI)
- `pygtgraph` (for plotting in GUI)
