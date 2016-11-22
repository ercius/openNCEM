'''
Tests for the eva.ring_diff module.
'''

import unittest
import numpy as np
import matplotlib.pyplot as plt
import os

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.math
import ncempy.eval.ring_diff


class test_ringdiff(unittest.TestCase):
    '''
    Test the evaluation of ring diffraction patterns.
    '''
    
    def test_eva(self):
        
        plt.close('all')        
        show=True
        
        
        ### test settings writing and reading
        if os.path.isfile('ncempy/test/resources/output/test_settings.emd'):
            os.remove('ncempy/test/resources/output/test_settings.emd')
            
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/test_settings.emd')
        
        ## minimum settings necessary, with many nones
        settings_minimum = { 'lmax_r': 10,
                 'lmax_thresh': 600,
                 'lmax_cinit': (984, 1032),
                 'lmax_range': (6e9, 8e9),
                 'plt_imgminmax': None,
                 'ns': (2,3,4),
                 'rad_rmax': None,
                 'rad_dr': None,
                 'rad_sigma': None,
                 'mask': None,
                 'fit_rrange': (1.5e9, 9.5e9),
                 'back_xs': (1.3e9, 1.5e9, 9.05e9),
                 'back_xswidth': 0.05e9,
                 'back_init': (1, 1.0e12, -1.0),
                 'fit_funcs': ('voigt',),
                 'fit_init': ( 5e10, 7.3e9, 1.1e7, 2.5e7 ),
                 'fit_maxfev': None
               }
            
        # put       
        par = femd.file_hdl.create_group('put_minimum')
        
        # wrong input
        with self.assertRaises(TypeError):
            ncempy.eval.ring_diff.put_settings( 'not_par', settings_minimum )
        with self.assertRaises(TypeError):
            ncempy.eval.ring_diff.put_settings( par, 42 )
              
        ncempy.eval.ring_diff.put_settings( par, settings_minimum )
        
        # get
        
        # wrong input
        with self.assertRaises(TypeError):
            no_settings = ncempy.eval.ring_diff.get_settings('not a group')        
        
        in_settings_minimum = ncempy.eval.ring_diff.get_settings( par['settings_ringdiffraction'] )
        
        # compare the two settings
        for key in settings_minimum:
            self.assertTrue(key in in_settings_minimum)
            
            comp = (settings_minimum[key]==in_settings_minimum[key])
            if isinstance(comp, bool):
                self.assertTrue(comp)
            else:
                self.assertTrue(comp.all())
        
        ## full settings with all parameters set
        settings_full = { 'lmax_r': 10,
                 'lmax_thresh': 600,
                 'lmax_cinit': (984, 1032),
                 'lmax_range': (6e9, 8e9),
                 'plt_imgminmax': (0.0, 0.2),
                 'ns': (2,3,4),
                 'rad_rmax': 42.,
                 'rad_dr': 0.1,
                 'rad_sigma': 0.01,
                 'mask': np.ones((25,25)),
                 'fit_rrange': (1.5e9, 9.5e9),
                 'back_xs': (1.3e9, 1.5e9, 9.05e9),
                 'back_xswidth': 0.05e9,
                 'back_init': (1, 1.0e12, -1.0),
                 'fit_funcs': ('voigt',),
                 'fit_init': ( 5e10, 7.3e9, 1.1e7, 2.5e7 ),
                 'fit_maxfev': 1000
               }
               
        # put       
        par = femd.file_hdl.create_group('put_full')
        ncempy.eval.ring_diff.put_settings( par, settings_full )
        
        # get
        in_settings_full = ncempy.eval.ring_diff.get_settings( par['settings_ringdiffraction'] )
        
        # compare the two settings
        for key in settings_full:
            self.assertTrue(key in in_settings_full)
            
            comp = (settings_full[key]==in_settings_full[key])
            if isinstance(comp, bool):
                self.assertTrue(comp)
            else:
                self.assertTrue(comp.all())
        
        
        ### single evaluation
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Pt_SAED_D910mm_single/Pt_SAED_D910mm_single.emd')
        emdgrp = femd.list_emds[0]
        
        if os.path.isfile('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_single.emd'):
            os.remove('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_single.emd')
        femd_out = ncempy.io.emd.fileEMD('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_single.emd')
        
        # write evaluation details
        grp_eva = femd_out.file_hdl.create_group('evaluation')
        
        # wrong input
        with self.assertRaises(TypeError):
            nohdl = ncempy.eval.ring_diff.put_sglgroup('nogroup', 'not_working', emdgrp)
        with self.assertRaises(TypeError):
            nohdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'not_working', 'nogroup')
        
        hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Pt_SAED_D910mm_single', emdgrp)

        # put the settings
        settings = { 'lmax_r': 10,
                 'lmax_thresh': 600,
                 'lmax_cinit': (984, 1032),
                 'lmax_range': (6e9, 8e9),
                 'plt_imgminmax': (0.,0.2),
                 'ns': (2,3,4),
                 'rad_rmax': None,
                 'rad_dr': None,
                 'rad_sigma': None,
                 'mask': None,
                 'fit_rrange': (1.5e9, 9.5e9),
                 'back_xs': (1.3e9, 1.5e9, 9.05e9),
                 'back_xswidth': 0.05e9,
                 'back_init': (1, 1.0e12, -1.0),
                 'fit_funcs': ('voigt',),
                 'fit_init': ( 5e10, 7.3e9, 1.1e7, 2.5e7 ),
                 'fit_maxfev': None
               }
        ncempy.eval.ring_diff.put_settings(femd_out.file_hdl, settings)

        # run the evaluation
        # wrong input
        with self.assertRaises(TypeError):
            ncempy.eval.ring_diff.run_sglgroup(42, femd_out)
        with self.assertRaises(TypeError):
            ncempy.eval.ring_diff.run_sglgroup(hdl, 'notanemdfile')
                    
        ncempy.eval.ring_diff.run_sglgroup(hdl, femd_out, verbose=True, showplots=False)
        
        
        ### evaluation of a series
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/Au_SAED_D910mm_20x_at_800.emd')
        emdgrp = femd.list_emds[0]
        
        if os.path.isfile('ncempy/test/resources/output/evaluation_Au_SAED_D910mm_20x_at_800.emd'):
            os.remove('ncempy/test/resources/output/evaluation_Au_SAED_D910mm_20x_at_800.emd')
        femd_out = ncempy.io.emd.fileEMD('ncempy/test/resources/output/evaluation_Au_SAED_D910mm_20x_at_800.emd')
        
        # write evaluation details
        grp_eva = femd_out.file_hdl.create_group('evaluation')
        
        hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Au_SAED_D910mm_20x_at_800', emdgrp)
        
        # put the settings
        settings = { 'lmax_r': 16,
                     'lmax_thresh': 250,
                     'lmax_cinit': (1012, 1020),
                     'lmax_range': (6.5e9, 7.3e9),
                     'plt_imgminmax': (0.,0.4),
                     'ns': (2,3,4),
                     'rad_rmax': None,
                     'rad_dr': None,
                     'rad_sigma': None,
                     'mask': None,
                     'fit_rrange': (1.5e9, 9.5e9),
                     'back_xs': (1.5e9, 5.7e9, 9.05e9),
                     'back_xswidth': 0.05e9,
                     'back_init': (1, 1.5e9, -0.8),
                     'fit_funcs': ('voigt','voigt','voigt','voigt','voigt'),
                     'fit_init': ( 3.4e10, 4.25e9, 1.1e7, 2.5e7,
                                   9.1e9, 4.90e9, 1.1e7, 2.5e7,
                                   1.6e10, 6.95e9, 1.1e7, 2.5e7,
                                   1.1e10, 8.14e9, 1.1e7, 2.5e7,
                                   1.6e9, 8.50e9, 1.1e7, 2.5e7 ),
                     'fit_maxfev': None
                   }
        ncempy.eval.ring_diff.put_settings(femd_out.file_hdl, settings)

        # run the evaluation
        ncempy.eval.ring_diff.run_sglgroup(hdl, femd_out, verbose=True, showplots=False)
        
        
        ### run multiple evaluations
        ### single evaluation
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Pt_SAED_D910mm_single/Pt_SAED_D910mm_single.emd')
        emdgrp = femd.list_emds[0]
        
        if os.path.isfile('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_multiple.emd'):
            os.remove('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_multiple.emd')
        femd_out = ncempy.io.emd.fileEMD('ncempy/test/resources/output/evaluation_Pt_SAED_D910mm_multiple.emd')
        
        # write evaluation details
        grp_eva = femd_out.file_hdl.create_group('evaluation')
        hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Pt_SAED_D910mm_first', emdgrp)
        hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Pt_SAED_D910mm_second', emdgrp)
        hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, 'Pt_SAED_D910mm_third', emdgrp)

        # put the settings
        settings = { 'lmax_r': 10,
                 'lmax_thresh': 600,
                 'lmax_cinit': (984, 1032),
                 'lmax_range': (6e9, 8e9),
                 'plt_imgminmax': (0.,0.2),
                 'ns': (2,3,4),
                 'rad_rmax': None,
                 'rad_dr': None,
                 'rad_sigma': None,
                 'mask': None,
                 'fit_rrange': (1.5e9, 9.5e9),
                 'back_xs': (1.3e9, 1.5e9, 9.05e9),
                 'back_xswidth': 0.05e9,
                 'back_init': (1, 1.0e12, -1.0),
                 'fit_funcs': ('voigt',),
                 'fit_init': ( 5e10, 7.3e9, 1.1e7, 2.5e7 ),
                 'fit_maxfev': None
               }
        ncempy.eval.ring_diff.put_settings(grp_eva, settings)
        
        ncempy.eval.ring_diff.run_all(grp_eva, femd_out, verbose=True, showplots=False)
        

        # wait for plots
        if show:
            plt.show()     
        

# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
