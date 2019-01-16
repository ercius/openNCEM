ncempy.io
==========

The ``ncempy.io`` module contains the file IO necessary for the various dataformats floating around in electron microscopy. Internally the package is designed to work in the versatile EMD file format, other formats are interfaced with importers. The interfaces of each file format are implemented in their own classes.

Contents
--------

Overview of contents with short description:

+--------------------+--------------------------------------------------------------------+
| Module             | Description                                                        |
+====================+====================================================================+
| emd                | EMD file format.                                                   |
+--------------------+--------------------------------------------------------------------+
| mrc                | MRC file format.                                                   |
+--------------------+--------------------------------------------------------------------+
| dm                 | DM3 and DM4 file format used by Gatan DigitalMicrograph.           |
+--------------------+--------------------------------------------------------------------+
| ser                | SER file format used by TIA (FEI/Thermo Fischer).                  |
+--------------------+--------------------------------------------------------------------+
| emdVelox           | HDF5 file format used by Velox (FEI/Thermo Fischer).               |
+--------------------+--------------------------------------------------------------------+
