"""
Tests for the algo.peak_find module.
"""

import pytest
from pathlib import Path

import numpy as np

import ncempy.io.emd
import ncempy.algo.peak_find


class Test_peak_find:
    """ Test the algo.peak_find module.

    """

    @pytest.fixture
    def data_location(self):
        """Get the location of the test data files"""
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    @pytest.fixture
    def one_peak_2d(self):
        """Generate an image with a single peak in it for simple testing"""
        im = np.zeros((64, 32), dtype=np.uint8)
        im[31, 15] = 1
        return im

    @pytest.fixture
    def one_peak_3d(self):
        """Generate an image with a single peak in it for simple testing"""
        im = np.zeros((64, 32, 16), dtype=np.uint8)
        im[31, 15, 7] = 1
        return im

    @pytest.fixture
    def one_gaussian_2d(self):
        """Generate an image with a single gaussian in it for simple testing"""
        peaks = np.asarray(((10, 11), (33, 14)))
        im = ncempy.algo.peak_find.peaksToImage(peaks, (64, 32), (1, 2), (5, 5))
        return im
    
    def test_doubleRoll(self, one_peak_2d):
        assert one_peak_2d[31, 15] == 1
        one_peak_roll = ncempy.algo.peak_find.doubleRoll(one_peak_2d, (1, 1))
        assert one_peak_roll[32, 16] == 1

    def test_tripleRoll(self, one_peak_3d):
        assert one_peak_3d[31, 15, 7] == 1
        one_peak_roll = ncempy.algo.peak_find.tripleRoll(one_peak_3d, (1, 1, 1))
        assert one_peak_roll[32, 16, 8] == 1

    def test_peakFind2D(self, one_gaussian_2d):
        peaks = ncempy.algo.peak_find.peakFind2D(one_gaussian_2d, 0.5)
        assert np.all(peaks == ((10, 11), (33, 14)))

    def test_lattice2D_norm(self):
        u = (1, 0)
        v = (0, 1)
        a = 1.5
        b = 1.75
        n = (3, 3)
        lat = ncempy.algo.peak_find.lattice2D_norm(u, v, a, b, (0, 0), n)
        assert np.allclose(lat[1, :], (1.5, 0))

