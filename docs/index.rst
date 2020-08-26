openNCEM Documentation
======================

A collection of packages and tools for electron microscopy data analysis supported by the National Center for Electron Microscopy facility of the Molecular Foundry.

Components
----------

The openNCEM collection comes with different components. The general functionality is condensed in libraries/packages, distinguished by programming language or framework. Tools are provided, addressing particular problems or tasks, leveraging the functionality provided in the libraries/packages.

.. toctree::
   :maxdepth: 2

    ncempy - openNCEM's python package<ncempy>
    tools - wrapping the provided functionality into useful tools<tools>


Examples
--------
Simple file loading
^^^^^^^^^^^^^^^^^^^
For simple file loading use the convenient Reader functions. An example of reading a DM file is below.

DM files
^^^^^^^^
Here is an example of loading a Digital Micrograph file. It returns the data and the metadata in
an easily human readable form. 

.. code-block:: python
   
   >>> from ncempy.io import dm; import matplotlib.pyplot as plt 
   >>> dmData = dm.dmReader('path/to/file/data.dm3') #a simple one image data file
   >>> plt.imshow(dmData['data']) #show the image using pyplot
   
For developers
^^^^^^^^^^^^^^
Developers will want access to the entire internal structure of the file.
Then you need to instantiate the class, parse the header and then access its contents

Again, here is a more thorough way to open a DM file
.. code-block:: python

   >>> from ncempy.io import dm; import matplotlib.pyplot as plt 
   >>> dm1 = dm.fileDM('/path/to/file/data.dm3')
   >>> dmData = dm1.getDataset(0) #the full first data set
   >>> dmSlice = dm1.getSlice(0,1) #the second image of the first data set; equal to dmData[1,:,:] from the line above
   >>> print(dm1.metaData) #all of the interesting metadata for this data set (pixel size, accelerating voltage, etc.)
   >>> print(dm1.allTags) #all of the tags in the file including tags specific to the Digital Micrpograph software program
    
    
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
