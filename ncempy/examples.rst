---------------
ncempy examples
---------------

Directory for examples within the NCEMpy project.

Jupyter notebooks
-----------------

* notebooks/3D_slicer.ipynb
    Uses ipywidgets to load a 3D dataset and allows the user to page through the different images.

* notebooks/ncempy io examples.ipynb
    Examples of loading data sets from DM3/4, SER, MRC and EMD files using the low- and high-level readers of the ncempy.io module.

* Velox EMD review
    An example of how to read EMD Velox images and browse their metadata.

EDSTomo
-------

* data/L2083-K-4-1/ConvertBrukerToEMD/ConvertBrukerToEMD.ipynb
    Example to convert a tomographic tilt series of Bruker bcf files into the NCEMpy EMD format.  This example has a dependency on the hyperspy package.


* data/L2083-K-4-1/ReconstructEDSTomo/ReconstructEDSTomo.ipynb
    Example to reconstruct an EDS tomographic tilt series using the EDSTomo module and GENFIRE.
