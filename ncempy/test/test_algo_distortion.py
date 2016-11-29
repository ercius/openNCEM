'''
Tests for the algo.distortion module.
'''

import unittest
import numpy as np
import matplotlib.pyplot as plt
import os

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.distortion

class test_ringdiff(unittest.TestCase):
    '''
    Test the distortion module on ring diffraction patterns.
    '''
   
    
    def test_distortion(self):
        '''
        Test the distortion fitting algorithms to be used on ring diffraction patterns.
        '''
        
        plt.close('all')
        show=True
        
        # get an image
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Pt_SAED_D910mm_single/Pt_SAED_D910mm_single.emd')
        data, dims = femd.get_emdgroup(femd.list_emds[0])
        
        img = data[0,:,:]
        dims = dims[1:3]
        
        # and a list of maxima
        points = ncempy.algo.local_max.local_max(img, 10, 600)
        points = ncempy.algo.local_max.points_todim(points, dims)
        
        # parameters
        center_init = ncempy.algo.local_max.points_todim((984,1032), dims)
        
        
        ## filter_ring
        # wrong input
        with self.assertRaises(TypeError):
            nopoints = ncempy.algo.distortion.filter_ring(np.array([[0,1,2,3,4],[5,6,7,8,9]]), center_init, (6e9, 8e9))
        with self.assertRaises(TypeError):
            nopoints = ncempy.algo.distortion.filter_ring(42, center_init, (6e9, 8e9))    
        with self.assertRaises(TypeError):
            nopoints = ncempy.algo.distortion.filter_ring(points, (1,2,3,4), (6e9, 8e9))        
            
        # no points
        nopoints = ncempy.algo.distortion.filter_ring(points, center_init, (8e9, 6e9))
        self.assertIsNone(nopoints)
        
        # feed through center
        self.assertTrue(np.array_equal(center_init, ncempy.algo.distortion.filter_ring(center_init, center_init, (0,0))))
        
        # working
        points = ncempy.algo.distortion.filter_ring(points, (center_init[0,0], center_init[0,1]), (6e9, 8e9))
        self.assertIsNotNone(points)
        points = ncempy.algo.distortion.filter_ring(points, center_init, (6e9, 8e9))
        self.assertIsNotNone(points)
        
        
        ## points_topolar
        # wrong input
        with self.assertRaises(TypeError):
            nopoints_plr = ncempy.algo.distortion.points_topolar(42, center_init)
        with self.assertRaises(TypeError):
            nopoints_plr = ncempy.algo.distortion.points_topolar(points, (1,2,3))
        with self.assertRaises(TypeError):
            nopoints_plr = ncempy.algo.distortion.points_topolar(np.array([[0,1,2],[3,4,5],[6,7,8]]), center_init)
        
        # working
        points_plr = ncempy.algo.distortion.points_topolar(points, center_init)
        
        
        ## plot_ringpolar
        # wrong input
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_ringpolar(np.array([0,1,2,3,4]), dims)
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_ringpolar(42, dims)
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_ringpolar(points_plr, dims[0])
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_ringpolar(points_plr, 42)
                    
        # working
        plot = ncempy.algo.distortion.plot_ringpolar(points_plr, dims, show=show)
        
        
        ## optimize_center
        # wrong input
        with self.assertRaises(TypeError):
            nocenter = ncempy.algo.distortion.optimize_center(np.array([0,1,2,3,4]), center_init)
        with self.assertRaises(TypeError):
            nocenter = ncempy.algo.distortion.optimize_center(42, center_init)
        with self.assertRaises(TypeError):
            nocenter = ncempy.algo.distortion.optimize_center(points, (1,2,3,4,5))
        with self.assertRaises(TypeError):
            nocenter = ncempy.algo.distortion.optimize_center(points, 'string')
        
        # get warning
        #center = ncempy.algo.distortion.optimize_center(points, center_init, maxfev=1)
        # working
        center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=True)
        
        points_plr = ncempy.algo.distortion.points_topolar(points, center)
        plot = ncempy.algo.distortion.plot_ringpolar(points_plr, dims, show=show)
    
        ## optimize_distortion
        # wrong input
        with self.assertRaises(TypeError):
            nodists = ncempy.algo.distortion.optimize_distortion(np.array([[0,1,2,3,4],[5,6,7,8]]), (2,3,4))
        with self.assertRaises(TypeError):
            nodists = ncempy.algo.distortion.optimize_distortion(42, (2,3,4))        
        with self.assertRaises(TypeError):
            nodists = ncempy.algo.distortion.optimize_distortion(points_plr, 2)
            
        # working
        dists = ncempy.algo.distortion.optimize_distortion(points_plr, (2,3,4))
        
        ## plot_distpolar
        # wrong input
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_distpolar(np.array([[0,1,2,3],[4,5,6,7]]), dims, dists, (2,3,4))
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_distpolar(42, dims, dists, (2,3,4))
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_distpolar(points_plr, dims[0], dists, (2,3,4))
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_distpolar(points_plr, dims, dists[0:3], (2,3,4))
        with self.assertRaises(TypeError):
            noplot = ncempy.algo.distortion.plot_distpolar(points_plr, dims, dists, 42)
            
        
        # try a bunch of fits
        for ns in ( (2,), (2,3), (2,3,4), (2,3,4,5), (4,) ):
     
            dists = ncempy.algo.distortion.optimize_distortion(points_plr, ns, verbose=True)
            plot = ncempy.algo.distortion.plot_distpolar(points_plr, dims, dists, ns, show=show)
            

        # wait for plots
        if show:
            plt.show()
    

# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
