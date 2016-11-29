'''
Tests for the algo.radial_profile module.
'''

import unittest
import numpy as np
import matplotlib.pyplot as plt
import os

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.distortion
import ncempy.algo.math
import ncempy.algo.radial_profile


class test_ringdiff(unittest.TestCase):
    '''
    Test the algo.radial_profile module.
    '''
    
    def test_radialprofile(self):
        '''
        Test the radial_profile algorithms to be used on ring diffraction patterns.
        '''
        
        plt.close('all')        
        show=True
        
        # get an image
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Pt_SAED_D910mm_single/Pt_SAED_D910mm_single.emd')
        data, dims = femd.get_emdgroup(femd.list_emds[0])
        
        img = data[0,:,:]
        dims = dims[1:3]
        
        # points on a ring
        points = ncempy.algo.local_max.local_max(img, 10, 600)
        points = ncempy.algo.local_max.points_todim(points, dims)
        center_init = ncempy.algo.local_max.points_todim((984,1032), dims)
        points = ncempy.algo.distortion.filter_ring(points, center_init, (6e9, 8e9))
        
        # optimize center
        center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=True)
        
        # fit distortions
        ns = (2,3,4)
        points_plr = ncempy.algo.distortion.points_topolar(points, center)
        dists = ncempy.algo.distortion.optimize_distortion(points_plr, ns)
        
        # check input
        plot = ncempy.algo.distortion.plot_distpolar(points_plr, dims, dists, ns, show=show)
        
        
        ## corrected images
        img_corr = ncempy.algo.radial_profile.correct_distortion( img, dims, center, ns, dists )
        
        
        ## calc_polarcoords
        # wrong input
        with self.assertRaises(TypeError):
            nors, nothes = ncempy.algo.radial_profile.calc_polarcoords( (1,2,3,4,5) )
        with self.assertRaises(TypeError):
            nors, nothes = ncempy.algo.radial_profile.calc_polarcoords( center, dims[0] )
        with self.assertRaises(TypeError):
            nors, nothes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, ns )            
        with self.assertRaises(TypeError):
            nors, nothes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, 42, dists )                    
        with self.assertRaises(TypeError):
            nors, nothes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, ns, dists[0:5] )

        # working
        rs_nodist, thes_nodist = ncempy.algo.radial_profile.calc_polarcoords( center, dims )
        rs, thes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, ns, dists )


        ## calc_radialprofile
        rMax = np.abs(dims[0][0][0]-dims[0][0][1])*np.min(img.shape)/2.0
        dr = np.abs(dims[0][0][0]-dims[0][0][1])/10.
        rsigma = np.abs(dims[0][0][0]-dims[0][0][1])
        # wrong input
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( 42, rs, rMax, dr, rsigma )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, 42, rMax, dr, rsigma )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs[0:25,0:25], rMax, dr, rsigma )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, 'one', dr, rsigma )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, 'two', rsigma )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, 'three' )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma, mask=42 )
        with self.assertRaises(TypeError):
            noR, noI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma, mask=np.ones((25,25)) )
                                                
        # working
        R, I = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma )
        
        # everything masked
        emR, emI = ncempy.algo.radial_profile.calc_radialprofile( img, rs, rMax, dr, rsigma, mask=np.zeros(img.shape) )


        ## plot_radialprofile
        # wrong input
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_radialprofile( 42, I, dims )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_radialprofile( R, 42, dims )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_radialprofile( R, I, 'notdims' )
        
        # working
        plot = ncempy.algo.radial_profile.plot_radialprofile( R, I, dims, show=show) 
        plot = ncempy.algo.radial_profile.plot_radialprofile( emR, emI, dims, show=show)
        
        
        ## compare without distortion correction
        R_nodist, I_nodist = ncempy.algo.radial_profile.calc_radialprofile( img, rs_nodist, rMax, dr, rsigma )
        plot = ncempy.algo.radial_profile.plot_radialprofile( R_nodist, I_nodist, dims, show=show)
        
        
        # cut radial profile
        sel = (R>=1.5e9)*(R<=9.5e9)
        I = I[sel]
        R = R[sel]
        

        ## fit_radialprofile
        funcs = [ 'const', 'powlaw', 'voigt' ]
        init_guess = (  10,
                        1.0e12, -1.0,
                        5e10, 7.3e9, 1.1e7, 2.5e7 )
        
        # wrong input
        with self.assertRaises(TypeError):
            nores = ncempy.algo.radial_profile.fit_radialprofile( 42, I, funcs, init_guess )
        with self.assertRaises(TypeError):
            nores = ncempy.algo.radial_profile.fit_radialprofile( R, 42, funcs, init_guess )
        with self.assertRaises(TypeError):
            nores = ncempy.algo.radial_profile.fit_radialprofile( R, I[0:25], funcs, init_guess )
        with self.assertRaises(TypeError):
            nores = ncempy.algo.radial_profile.fit_radialprofile( R, I, funcs, init_guess[0:3] )
        with self.assertRaises(TypeError):
            nores = ncempy.algo.radial_profile.fit_radialprofile( R, I, 42, init_guess )
        
        # working
        res = ncempy.algo.radial_profile.fit_radialprofile( R, I, funcs, init_guess, maxfev=10000 ) 
                 
        
        ## plot_fit
        # wrong input
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_fit( 42, I, dims, funcs, res )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_fit( R, 42, dims, funcs, res )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_fit( R, I, 'notdims', funcs, res )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs, res[0:3] )
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.radial_profile.plot_fit( R, I, dims, 42, res )

        # before fit
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs, init_guess, show=show )
        # after fit
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs, res, show=show )
        
        
        # try two voigts
        funcs = [ 'voigt', 'voigt' ]
        init_guess = ( 6e10, 7.3e9, 2e7, 2e7, 
                       8e10, 2.4e9, 1e8, 1e8 )
        res = ncempy.algo.radial_profile.fit_radialprofile( R, I, funcs, init_guess, maxfev=10000 )
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs, res, show=show )


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
            
        funcs = [ 'const', 'powlaw']
        init_guess = (1, 1.0e12, -1.0)
        res = ncempy.algo.radial_profile.fit_radialprofile( fit_R, fit_I, funcs, init_guess, maxfev=1000 )
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs, res, show=show )
        
        I = I - ncempy.algo.math.sum_functions(R, funcs, res)
        
        plot = ncempy.algo.radial_profile.plot_radialprofile( R, I, dims, show=show )

        #import pdb;pdb.set_trace()
        
        # wait for plots
        if show:
            plt.show()            
    

# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
