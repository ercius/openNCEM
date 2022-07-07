ncempy.algo
===========

The ``ncempy.algo`` module contains the algorithms used for the various image processing, data evaluating and simulating tasks. These are to be reused for different evaluations and are therefore written in a general way.
The evaluation targeting a certain problem is supposed to go in ``ncempy.eval`` and import all the necessary algorithms from ``ncempy.algo``.

Contents
--------

Overview of contents with short description:

+--------------------+--------------------------------------------------------------------+
| Module             | Description                                                        |
+====================+====================================================================+
| distortion         | Treat distortion in diffraction patterns.                          |
+--------------------+--------------------------------------------------------------------+
| local_max          | Find local maxima in an image.                                     |
+--------------------+--------------------------------------------------------------------+
| math               | Flexible fit function construction.                                |
+--------------------+--------------------------------------------------------------------+
| radial_profile     | Calculate azimuthally averaged radial profiles of images.          |
+--------------------+--------------------------------------------------------------------+
| multicorr_funcs    | Helper functions for multicorr phase and cross correlation         |
+--------------------+--------------------------------------------------------------------+
| peak_find          | Find peaks in images, fit to a lattice, calcualte unit calls       |
+--------------------+--------------------------------------------------------------------+
| gaussND            | Multi-dimensional Gaussian (and Lorentz) functions for fitting.    |
+--------------------+--------------------------------------------------------------------+
