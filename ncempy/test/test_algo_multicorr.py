'''
Tests for the algo.multicorr module.
'''

import unittest
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from io import StringIO
import cv2

import ncempy.algo.multicorr as mc

class test_shiftedimages(unittest.TestCase):
    """Test the multicorr module on two shifted images"""

    def test_parse_input(self):
        '''
        Test to make sure inputs are properly accounted for.
        '''


        g1 = 1
        g2 = np.zeros([4,4])
        with self.assertRaises(TypeError):
            mc.multicorr(g1, g2).parse_input() # both must be arrays
        with self.assertRaises(TypeError):
            mc.multicorr(g2, g1).parse_input()

        # tests input sanitation.
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out
            mc.multicorr(g2, g2, 'dummy', 4).parse_input()
            output = out.getvalue().strip()
            self.assertEqual(output, 'Unknown method used, setting to cross')

            out = StringIO()
            sys.stdout = out
            mc.multicorr(g2, g2, 'cross', 23.1).parse_input()
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer, rounding down')

            out = StringIO()
            sys.stdout = out
            mc.multicorr(g2, g2, 'cross', -23.1).parse_input()
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer, rounding down\nUpsample factor is < 1, setting to 1')

            out = StringIO()
            sys.stdout = out
            mc.multicorr(g2, g2, 'cross', 'phase').parse_input()
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer or float, setting to 1')
        finally:
            sys.stdout = saved_stdout

        with self.assertRaises(TypeError):
            mc.multicorr(np.zeros([2,3]), np.zeros([3,2]), 'phase', 1).parse_input()

    def test_multicorr_math_initial_correlation_image(self):
        '''
        Test to check if the correlation is working.
        '''
        g2 = np.zeros((3,3))
        test = mc.multicorr(g2, g2)
        assert test.G2.all() == g2.all()

        filename = '/Users/Tom/Downloads/matt_beard.jpg'
        filename_shifted = '/Users/Tom/Downloads/matt_beard_shifted.jpg'
        G1 = cv2.imread(filename, 0)
        G2 = cv2.imread(filename_shifted, 0)
        out_phase = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'phase', 2)
        out_phase.multicorr()

        np.testing.assert_almost_equal(np.exp(1j * np.angle((np.multiply(np.fft.fft2(G1), np.conj(np.fft.fft2(G2)))))), out_phase.imageCorr, decimal=4)

        np.testing.assert_almost_equal(np.real(np.fft.ifft2(out_phase.imageCorr)), out_phase.imageCorrIFT, decimal=4)

        # plt.imshow(out.imageCorrIFT)
        # plt.show(block = True)
        self.assertEqual(out_phase.xyShift, [-30.0, 0.0])

        # test cross correlation
        out_cross = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 1)
        out_cross.multicorr()

        np.testing.assert_almost_equal(np.multiply(np.fft.fft2(G1), np.conj(np.fft.fft2(G2))), out_cross.imageCorr, decimal=4)

        # test hybrid correlation and dft upsample
        out_hybrid = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
        out_hybrid.multicorr()
        print(type(out_hybrid.xyShift))
        self.assertEqual(list(out_hybrid.xyShift), [-30.0, 0.0])

        # plt.imshow(np.subtract(G1, np.real(np.fft.ifft2(out.G2shift))))
        # plt.show(block = True)



if __name__ == '__main__':
    unittest.main()
