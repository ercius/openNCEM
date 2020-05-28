"""
Tests for the algo.radial_profile module.
"""

import pytest
from pathlib import Path

import numpy as np

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.distortion
import ncempy.algo.math
import ncempy.algo.radial_profile


class TestRingdiff():
    """
    Test the algo.radial_profile module.

    This mostly tests that things run without errors. Output
    tests have been removed.
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_radialprofile(self, data_location):
        """
        Test the radial_profile algorithms to be used on ring diffraction patterns.
        """

        with ncempy.io.emd.fileEMD(data_location / Path('Pt_SAED_D910mm_single.emd')) as emd0:
            data, dims = emd0.get_emdgroup(0)

            # Data comes out as 3D, but its really 1D
            img = data[0, :, :]
            dims = dims[1:3]
        
            # points on a ring
            points = ncempy.algo.local_max.local_max(img, 10, 600)
            points = ncempy.algo.local_max.points_todim(points, dims)
            center_init = ncempy.algo.local_max.points_todim(np.array((984, 1032)), dims)
            points = ncempy.algo.distortion.filter_ring(points, center_init, (6e9, 8e9))

            # optimize center
            center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=True)

            # fit distortions
            ns = (2, 3, 4)
            points_plr = ncempy.algo.distortion.points_topolar(points, center)
            dists = ncempy.algo.distortion.optimize_distortion(points_plr, ns)

            # corrected images
            img_corr = ncempy.algo.radial_profile.correct_distortion( img, dims, center, ns, dists )

            # working
            rs_nodist, thes_nodist = ncempy.algo.radial_profile.calc_polarcoords( center, dims )
            rs, thes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, ns, dists )

            # calc_radialprofile
            rMax = np.abs(dims[0][0][0]-dims[0][0][1])*np.min(img.shape)/2.0
            dr = np.abs(dims[0][0][0]-dims[0][0][1])/10.
            rsigma = np.abs(dims[0][0][0]-dims[0][0][1])

            # working
            R, I = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma )

            # everything masked
            emR, emI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma, mask=np.zeros(img.shape) )

            ## compare without distortion correction
            R_nodist, I_nodist = ncempy.algo.radial_profile.calc_radialprofile(img, rs_nodist, rMax, dr, rsigma)

            # cut radial profile
            sel = (R >= 1.5e9) * (R <= 9.5e9)
            I = I[sel]
            R = R[sel]

            # fit_radialprofile
            funcs = ['const', 'powlaw', 'voigt']
            init_guess = (  10,
                            1.0e12, -1.0,
                            5e10, 7.3e9, 1.1e7, 2.5e7 )

            res = ncempy.algo.radial_profile.fit_radialprofile( R, I, funcs, init_guess, maxfev=10000 )

            # try two voigts
            funcs = ['voigt', 'voigt']
            init_guess = ( 6e10, 7.3e9, 2e7, 2e7,
                           8e10, 2.4e9, 1e8, 1e8 )
            res = ncempy.algo.radial_profile.fit_radialprofile(R, I, funcs, init_guess, maxfev=10000)

            # subtract a power law background fitted to specific points
            # get some fixpoints
            fit_xs = [1.3e9, 1.5e9, 9.05e9]
            fit_xswidth = 0.05e9
            fit_R = np.array([])
            fit_I = np.array([])
            for xpoint in fit_xs:
                ix = np.where( np.abs(R-xpoint) < fit_xswidth )
                fit_R = np.append(fit_R, R[ix])
                fit_I = np.append(fit_I, I[ix])

            funcs = ['const', 'powlaw']
            init_guess = (1, 1.0e12, -1.0)
            res = ncempy.algo.radial_profile.fit_radialprofile( fit_R, fit_I, funcs, init_guess, maxfev=1000 )

            I = I - ncempy.algo.math.sum_functions(R, funcs, res)
