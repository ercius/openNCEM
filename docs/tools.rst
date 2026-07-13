tools
=====

The collection of tools provided with openNCEM. These are based on the algorithms and routines contained in ``ncempy`` and are meant to provide solutions to specific problems. They are subordered by the experiment/simulation/task they address.

The GUI tools come with additional requirements:

* PySide (for GUI)
* pyqtgraph (for plotting in GUI)

RingDiff
--------

Tool to evaluate ring diffraction patterns. The command line tool is split into two parts: ``cmdline_ringdiff_prepare`` creates an evaluation file, and ``cmdline_ringdiff_run`` executes the evaluation. This allows settings to be modified with an external HDF5 tool (such as the HDF Viewer) before running the evaluation.

The GUI tool allows a standalone evaluation of ring diffraction patterns. It can read EMD files containing a single or a series of diffraction patterns, or read the settings and linked data from an evaluation file.

ser2emd
-------

A converter to read SER files and save them as EMD files. Optionally an EMI file is parsed as well to retrieve the metadata.

A command line as well as a GUI version are available. Run ``python3 gui_ser2emd.py`` for the GUI, or ``python3 cmdline_ser2emd.py --help`` for command line usage.

tif2emd
-------

GUI tools to convert TIF files to EMD files, for either a single TIF file or a series of TIF files.
