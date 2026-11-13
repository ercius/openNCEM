python_titan
============

A set of open source automation packages for Thermo Fisher (FEI) transmission electron microscopes. These scripts connect directly to the microscope and TIA using COM (via the ``comtypes`` package) to control acquisition and stage movements, and are not part of the installable ``ncempy`` package.

They currently support the NCEM TEAM Stage and the FEI Compustage.

STEM_minimal.py
----------------

A minimal script demonstrating how to connect to the microscope and TIA in STEM mode and acquire a single image (not saved).

STEMexperiments.py
-------------------

Acquires STEM data experiments such as focal series, scanning rotation series, time series, and single images. All data is saved as a Berkeley EMD file, readable with ``ncempy.io.emd``.

STEMTomo7.py
-------------

Acquires STEM tomography tilt series manually, saving all metadata and images to a Berkeley EMD file. Supports multiple images per tilt angle or a scan rotation series for dose fractionation and drift correction. Works with both the TEAM Stage and the FEI Compustage.

KStilter.py
------------

Determines the tilt required to bring a crystal on axis by clicking a position in a CBED pattern. Currently only works with the NCEM TEAM Stage and the FEI Flucam.
