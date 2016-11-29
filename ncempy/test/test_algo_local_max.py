'''
Tests for the algo.local_max module.
'''

import unittest
import numpy as np
import matplotlib.pyplot as plt
import os

import ncempy.io.emd
import ncempy.algo.local_max


class test_localmax(unittest.TestCase):
    '''
    Test the localmax module on diffraction patterns.
    '''
    
    def test_local_max(self):
        '''
        Test the local_maxima algorithm to be used for ring diffraction patterns.
        '''
        
        plt.close('all')        
        show=True
        
        # get an image
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Pt_SAED_D910mm_single/Pt_SAED_D910mm_single.emd')
        data, dims = femd.get_emdgroup(femd.list_emds[0])
        
        img = data[0,:,:]
        dims = dims[1:3]
        
        
        ## local_max
        # not working
        with self.assertRaises(TypeError):
            points = ncempy.algo.local_max.local_max(42, 10, 600)
        with self.assertRaises(TypeError):
            points = ncempy.algo.local_max.local_max(img, 'ten', 600)
        with self.assertRaises(TypeError):
            points = ncempy.algo.local_max.local_max(img, 10, 'sixhundred')
            
        # no points detected
        self.assertIsNone (ncempy.algo.local_max.local_max(img, 10, 100000))
        
        # working
        points = ncempy.algo.local_max.local_max(img, 10, 600)
        self.assertIsNotNone(points)
        
        
        ## plot_points
        # not working
        with self.assertRaises(TypeError):
            plot = ncempy.algo.local_max.plot_points(42, points)
        with self.assertRaises(TypeError):
            plot = ncempy.algo.local_max.plot_points(img, np.array([[0,1,2,3,4],[5,6,7,8,9]]))
            
        # px coords, non inverted
        plot = ncempy.algo.local_max.plot_points(img, points, vminmax=(0.0,0.2), show=show)
        
        
        ## points_todim
        # not working
        with self.assertRaises(TypeError):
            points = ncempy.algo.local_max.points_todim([0,1,2,3,5], dims)
        with self.assertRaises(TypeError):
            points = ncempy.algo.local_max.points_todim(points, dims[0])
        
        # convert points to real dim
        point = ncempy.algo.local_max.points_todim((0,0), dims)
        points = ncempy.algo.local_max.points_todim(points, dims)
        
        # real coords, inverted
        plot = ncempy.algo.local_max.plot_points(img, points, vminmax=(0.0,0.2), dims=dims, invert=True, show=show)
        
        
        # wait for plots
        if show:
            plt.show()
    
    
# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
