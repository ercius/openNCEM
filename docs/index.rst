openNCEM Documentation
======================

A collection of packages and tools for electron microscopy data analysis supported by the National Center for Electron Microscopy facility of the Molecular Foundry.

Installation
------------

We recommend installing the ``ncempy`` package from PyPi:

.. code-block:: bash

    >> pip install ncempy

If you wish to install the optional EDSTomo module then run:

.. code-block:: bash

    pip install 'ncempy[edstomo]'

Components
----------

The openNCEM collection comes with different components. The general functionality is condensed in libraries/packages, distinguished by programming language or framework. Tools are provided, addressing particular problems or tasks, leveraging the functionality provided in the libraries/packages.

.. toctree::
   :maxdepth: 2

    ncempy - openNCEM's python package<ncempy>
    tools - wraps the provided functionality in ncempy into useful tools<tools>
    python_titan - a set of open source automation packages for Thermo Fisher transmission electron microscopes<python_titan>

Examples
--------
Simple file reading
^^^^^^^^^^^^^^^^^^^
For simple file loading use the convenient `read` function.

Here is an example of loading a Digital Micrograph file. It returns the data and the metadata in
an easily human readable form. 

.. code-block:: python
   
   >>> import ncempy
   >>> data = ncempy.read('path/to/file/data.dm3')
   >>> print(data['data'].shape) # the shape of the data
   >>> print(data['pixelSize']) # print the pixel size

For developers
^^^^^^^^^^^^^^
Developers may want access to the entire internal structure of the file.
Instantiate the class, parse the header, and access the full contents.

Here is how to open a DM file and have access to the internal file parameters

.. code-block:: python

   >>> import ncempy
   >>> with ncempy.io.dm.fileDM('/path/to/file/data.dm3') as dm0:
   >>>     dmData = dm0.getDataset(0) #the full first data set
   >>>     dmSlice = dm0.getSlice(0, 1) #the second image of the first data set; equal to dmData[1,:,:] from the line above
   >>>     print(dm0.metaData) # all of the interesting metadata for this data set (pixel size, accelerating voltage, etc.)
   >>>     print(dm0.allTags) # all of the tags in the file including tags specific to the Digital Micrpograph software program

Reading Dectris Arina 4D-STEM data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Dectris Arina detector writes its 4D-STEM data set across an HDF5 master
file and one or more linked data files. Use `ncempy.io.dectris` to load the
full 4D data set from the master file.

.. code-block:: python

   >>> import ncempy
   >>> with ncempy.io.dectris.fileDECTRIS('/path/to/file/master.h5') as f0:
   >>>     data = f0.getDataset() # the full 4D data set, shape [scanY, scanX, frameY, frameX]
   >>>     print(data['data'].shape)
   >>>     print(f0.getMetadata()) # scan metadata (e.g. pixel size), if available

License
-------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.


    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


    You should have received a copy of the GNU General Public License
    along with this program.  If not, see https://www.gnu.org/licenses/.
