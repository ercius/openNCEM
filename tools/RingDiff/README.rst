--------
RingDiff
--------

Tool to evaluate ring diffraction patterns. The evaluation follows the routine described in F. Niekiel et al., submitted to Ultramicroscopy.

Command line usage
------------------

The command line tool is splitted into two parts. ``cmdline_ringdiff_prepare`` is used to create an evaluation file. This design allows to modify the various settings using an external HDF5 tool like the HDF Viewer, before executing the evaluation with the ``cmdline_ringdiff_run`` tool.

GUI tool
--------

The GUI tool allows a standalone evaluation of ring diffraction patterns. It can read EMD files containing a single or a series of diffraction patterns. Alternatively it can read the settings and the linked data from an evaluation file. To operate the evaluation follow the workflow on the left hand side from top to bottom and read what is happening in the log window.
