"""
Tests for the eva.ring_diff module.
"""

import pytest
from pathlib import Path
import tempfile

import numpy as np
import os

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.math
import ncempy.eval.ring_diff


class TestRingdiff():
    """
    Test the evaluation of ring diffraction patterns.
    """

    @pytest.fixture
    def temp_file(self):
        tt = tempfile.NamedTemporaryFile(mode='wb')
        tt.close()  # need to close the file to use it later
        return Path(tt.name)

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_eval_minimum(self, temp_file):
        with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:

            # minimum settings necessary, with many nones
            settings_minimum = {'lmax_r': 10,
                                'lmax_thresh': 600,
                                'lmax_cinit': (984, 1032),
                                'lmax_range': (6e9, 8e9),
                                'plt_imgminmax': None,
                                'ns': (2, 3, 4),
                                'rad_rmax': None,
                                'rad_dr': None,
                                'rad_sigma': None,
                                'mask': None,
                                'fit_rrange': (1.5e9, 9.5e9),
                                'back_xs': (1.3e9, 1.5e9, 9.05e9),
                                'back_xswidth': 0.05e9,
                                'back_init': (1, 1.0e12, -1.0),
                                'fit_funcs': ('voigt',),
                                'fit_init': (5e10, 7.3e9, 1.1e7, 2.5e7),
                                'fit_maxfev': None
                                }

            # put
            par = emd0.file_hdl.create_group('put_minimum')

            ncempy.eval.ring_diff.put_settings(par, settings_minimum)

            in_settings_minimum = ncempy.eval.ring_diff.get_settings(par['settings_ringdiffraction'])

            # compare the two settings
            for key in settings_minimum:
                assert (key in in_settings_minimum) is True

                comp = (settings_minimum[key] == in_settings_minimum[key])
                if isinstance(comp, bool):
                    assert comp is True
                #else:
                #    assert comp.all() is True

    def test_eval_full(self, temp_file):
        with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:
            # full settings with all parameters set
            settings_full = {'lmax_r': 10,
                             'lmax_thresh': 600,
                             'lmax_cinit': (984, 1032),
                             'lmax_range': (6e9, 8e9),
                             'plt_imgminmax': (0.0, 0.2),
                             'ns': (2, 3, 4),
                             'rad_rmax': 42.,
                             'rad_dr': 0.1,
                             'rad_sigma': 0.01,
                             'mask': np.ones((25, 25)),
                             'fit_rrange': (1.5e9, 9.5e9),
                             'back_xs': (1.3e9, 1.5e9, 9.05e9),
                             'back_xswidth': 0.05e9,
                             'back_init': (1, 1.0e12, -1.0),
                             'fit_funcs': ('voigt',),
                             'fit_init': (5e10, 7.3e9, 1.1e7, 2.5e7),
                             'fit_maxfev': 1000
                             }

            par = emd0.file_hdl.create_group('put_full')
            ncempy.eval.ring_diff.put_settings(par, settings_full)

            # get settings
            in_settings_full = ncempy.eval.ring_diff.get_settings(par['settings_ringdiffraction'])

            # compare the two settings
            for key in settings_full:
                assert (key in in_settings_full) is True

                comp = (settings_full[key] == in_settings_full[key])
                if isinstance(comp, bool):
                    assert comp is True
                #else:
                #    assert comp.all() is True

    def test_eval_single(self, data_location, temp_file):
        # single evaluation
        with ncempy.io.emd.fileEMD(data_location / Path('Pt_SAED_D910mm_single.emd')) as emd0:
            emdgrp = emd0.list_emds[0]

            with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd1:
                # write evaluation details
                grp_eva = emd1.file_hdl.create_group('evaluation')

                hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Pt_SAED_D910mm_single', emdgrp)

                # put the settings
                settings = {'lmax_r': 10,
                            'lmax_thresh': 600,
                            'lmax_cinit': (984, 1032),
                            'lmax_range': (6e9, 8e9),
                            'plt_imgminmax': (0., 0.2),
                            'ns': (2, 3, 4),
                            'rad_rmax': None,
                            'rad_dr': None,
                            'rad_sigma': None,
                            'mask': None,
                            'fit_rrange': (1.5e9, 9.5e9),
                            'back_xs': (1.3e9, 1.5e9, 9.05e9),
                            'back_xswidth': 0.05e9,
                            'back_init': (1, 1.0e12, -1.0),
                            'fit_funcs': ('voigt',),
                            'fit_init': (5e10, 7.3e9, 1.1e7, 2.5e7),
                            'fit_maxfev': None
                            }
                ncempy.eval.ring_diff.put_settings(emd1.file_hdl, settings)

                # run the evaluation
                ncempy.eval.ring_diff.run_sglgroup(hdl, emd1, verbose=True, showplots=False)
