ncempy.test
===========

This module contains the unittests for all other modules, functions, lines of code.

Requirements
------------

The tests are run on real data contained in a subfolder called ``resources``, the files have a size of several GBs. These are not included and the repo and have to be obtained elsewhere (for now just ask).

Usage
-----

For example:

``nosetests3 -s --with-coverage --cover-package=ncempy ncempy.test.test_io_emd``
