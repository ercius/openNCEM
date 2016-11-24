ncempy.algo
===========

The ``ncempy.algo`` module contains the algorithms used for the various image processing, data evaluating and simulating tasks. These are to be reused for different evaluations and are therefore written in a general way.
The evaluation targeting a certain problem is supposed to go in ``ncempy.eval`` and import all the necessary algorithms from ``ncempy.algo``.

+---------------------------------------+--------------------------------------------------------------------+
| Module                                | Description                                                        |
+=======================================+====================================================================+
| :doc:`ncempy.algo.distortion`         | Treat distortion in diffraction patterns.                          |
+---------------------------------------+--------------------------------------------------------------------+
| :doc:`ncempy.algo.local_max`          | Find local maxima in an image.                                     |
+---------------------------------------+--------------------------------------------------------------------+
| :doc:`ncempy.algo.math`               | Flexible fit function construction.                                |
+---------------------------------------+--------------------------------------------------------------------+
| :doc:`ncempy.algo.radial_profile`     | Calculate azimuthally averaged radial profiles of images.          |
+---------------------------------------+--------------------------------------------------------------------+

Contents
--------

.. toctree::
    :maxdepth: 1
    
    ncempy.algo.distortion
    ncempy.algo.local_max
    ncempy.algo.math
    ncempy.algo.radial_profile


