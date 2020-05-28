"""
Tests for the algo.distortion module.
"""

import pytest

from pathlib import Path

import numpy as np

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.distortion


class Testringdiff():
    """
    Test the distortion module on ring diffraction patterns.
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_distortion(self, data_location):
        """
        Test the distortion fitting algorithms to be used on ring diffraction patterns.
        """
        with ncempy.io.emd.fileEMD(data_location / Path('Pt_SAED_D910mm_single.emd')) as emd0:
            data, dims = emd0.get_emdgroup(0)

            # Data comes out as 3D, but its really 1D
            img = data[0, :, :]
            dims = dims[1:3]

            # and a list of maxima
            points = ncempy.algo.local_max.local_max(img, 10, 600)
            points = ncempy.algo.local_max.points_todim(points, dims)

            # parameters
            center_init = ncempy.algo.local_max.points_todim(np.array((984, 1032)), dims)

            # no points
            nopoints = ncempy.algo.distortion.filter_ring(points, center_init, (8e9, 6e9))
            assert nopoints is None

            # feed through center
            assert np.array_equal(center_init, ncempy.algo.distortion.filter_ring(center_init, center_init, (0, 0))) \
                   is True

            # working
            points = ncempy.algo.distortion.filter_ring(points, (center_init[0, 0], center_init[0, 1]), (6e9, 8e9))
            assert points is not None
            points = ncempy.algo.distortion.filter_ring(points, center_init, (6e9, 8e9))
            assert points is not None

            # points_topolar
            points_plr = ncempy.algo.distortion.points_topolar(points, center_init)
            assert points_plr is not None

            # optimize_center
            center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=True)

            points_plr = ncempy.algo.distortion.points_topolar(points, center)

            # optimize_distortion
            dists = ncempy.algo.distortion.optimize_distortion(points_plr, (2, 3, 4))

            # plot_distpolar
            # try a bunch of fits
            for ns in ((2,), (2, 3), (2, 3, 4), (2, 3, 4, 5), (4,)):
                dists = ncempy.algo.distortion.optimize_distortion(points_plr, ns, verbose=True)
