"""
Tests for the algo.local_max module.
"""

import pytest
from pathlib import Path

import numpy as np

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.viz


class TestLocalmax():
    """
    Test the localmax module on diffraction patterns.
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_local_max(self, data_location):
        """
        Test the local_maxima algorithm to be used for ring diffraction patterns.
        """
        
        with ncempy.io.emd.fileEMD(data_location / Path('Pt_SAED_D910mm_single.emd')) as emd0:
            data, dims = emd0.get_emdgroup(0)

            # Data comes out as 3D, but its really 1D
            img = data[0, :, :]
            dims = dims[1:3]

            # no points detected
            no_points = ncempy.algo.local_max.local_max(img, 10, 100000)
            assert no_points is None

            # working
            points = ncempy.algo.local_max.local_max(img, 10, 600)
            assert points is not None

            # convert points to real dim
            point = ncempy.algo.local_max.points_todim(np.array((0, 0)), dims)
            points = ncempy.algo.local_max.points_todim(points, dims)
