ncempy.io
==========

The ``ncempy.io`` module contains the file IO necessary for the various dataformats used in the electron microscopy
community.

The general ``read`` function will simply load any of the supported formats.

Each format can also be accessed using the low level interfaces of each file as a separate class.

Contents
--------

Overview of contents with short description:

+--------------------+--------------------------------------------------------------------+
| Module             | Description                                                        |
+====================+====================================================================+
| read               | Load any of the file formats below.                                |
+--------------------+--------------------------------------------------------------------+
| emd                | Berkeley EMD file format.                                          |
+--------------------+--------------------------------------------------------------------+
| mrc                | MRC file format which includes REC, ALI, and ST.                   |
+--------------------+--------------------------------------------------------------------+
| dm                 | DM3 and DM4 file format used by Gatan DigitalMicrograph.           |
+--------------------+--------------------------------------------------------------------+
| ser                | SER file format used by TIA (FEI/Thermo Fischer).                  |
+--------------------+--------------------------------------------------------------------+
| emdVelox           | HDF5 file format used by Velox (FEI/Thermo Fischer).               |
+--------------------+--------------------------------------------------------------------+
