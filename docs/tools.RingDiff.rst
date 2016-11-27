RingDiff
========

Tool to evaluate ring diffraction patterns. The evaluation follows the routine described in F. Niekiel et al., submitted to Ultramicroscopy.

Command line usage
------------------

The command line tool is splitted into two parts. ``cmdline_ringdiff_prepare`` is used to create an evaluation file. This design allows to modify the various settings using an external HDF5 tool like the HDF Viewer, before executing the evaluation with the ``cmdline_ringdiff_run`` tool.

GUI tool
--------

The GUI tool allows a standalone evaluation of ring diffraction patterns. It can read EMD files containing a single or a series of diffraction patterns. Alternatively it can read the settings and the linked data from an evaluation file. To operate the evaluation follow the workflow on the left hand side from top to bottom and read what is happening in the log window.

.. image:: https://cloud.githubusercontent.com/assets/23124266/20652218/dbff2964-b4a9-11e6-9985-b7e6594e2924.png
    :scale: 30 %

.. image:: https://cloud.githubusercontent.com/assets/23124266/20652219/de507d12-b4a9-11e6-8713-91211197aea9.png
    :scale: 30 %

.. image:: https://cloud.githubusercontent.com/assets/23124266/20652220/e02cd0c2-b4a9-11e6-9077-7d637fc2e19d.png
    :scale: 30 %

Requirements
^^^^^^^^^^^^

The GUI tool comes with additional requirements:

    * PySide (for GUI)
    * pygtgraph (for plotting in GUI)
