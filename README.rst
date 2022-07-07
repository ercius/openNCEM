========
openNCEM
========

A collection of packages, scripts, and tools for electron microscopy data analysis supported by the National Center for Electron Microscopy facility of the Molecular Foundry.

Read the documentation: https://openncem.readthedocs.io/en/latest/.

Installation of the ncempy package
==============================

The newest version can be installed using pip:

``>> pip install ncempy``

To ensure dependencies for the EDS tomography reconstruction are installed:

``>> pip install ncempy[edstomo]``

Components
==========

The openNCEM collection comes with different components, described below.

**ncempy**
    openNCEM's python package to provide various algorithms and routines to process or simulate images.


**tools**
    The tools build leveraging the algorithms and routines provided in the libraries/packages. They are subordered by the different problems they address.

Commands
========

**ncem2png**
    Extracts a PNG file from the data in a SER, DM3, or DM4 file. For 3D and 4D
    files, it tries to extract an image from the middle of the third and forth
    dimension.

License
=======

    ``ncempy`` is dual licensed under GPLv3 and MIT. The io module is the only part
    released under MIT to improve interoperability with other packages.

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
