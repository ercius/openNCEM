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
            mc.parse_input(g1, g2) # both must be arrays
        with self.assertRaises(TypeError):
            mc.parse_input(g2, g1)

        # tests input sanitation.
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out
            mc.parse_input(g2, g2, 'dummy', 4)
            output = out.getvalue().strip()
            self.assertEqual(output, 'Unknown method used, setting to cross')

            out = StringIO()
            sys.stdout = out
            mc.parse_input(g2, g2, 'cross', 23.1)
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer, rounding down')

            out = StringIO()
            sys.stdout = out
            mc.parse_input(g2, g2, 'cross', -23.1)
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer, rounding down\nUpsample factor is < 1, setting to 1')

            out = StringIO()
            sys.stdout = out
            mc.parse_input(g2, g2, 'cross', 'phase')
            output = out.getvalue().strip()
            self.assertEqual(output, 'Upsample factor is not an integer or float, setting to 1')
        finally:
            sys.stdout = saved_stdout

        with self.assertRaises(TypeError):
            mc.parse_input(np.zeros([2,3]), np.zeros([3,2]), 'phase', 1)

        with self.assertRaises(TypeError):
            mc.initial_correlation_image(g2, g2, 'tester', 1)

    def test_multicorr(self):
        '''
        Test to check if the correlation is working.
        '''
        g2 = np.zeros((3,3))
        test = mc.multicorr(np.fft.fft2(g2), np.fft.fft2(g2), 'phase', 3)
        print(test)
        assert list(test) == [0, 0]

        # filename = '/Users/Tom/Downloads/matt_beard.jpg'
        # filename_shifted = '/Users/Tom/Downloads/matt_beard_shifted.jpg'
        # G1 = cv2.imread(filename, 0)
        # G1 = G1[0:100, 0:100]
        G1 = np.zeros((100,100))
        G1[42,35] = 12
        # G2 = cv2.imread(filename_shifted, 0)
        G2 = np.real(np.fft.ifft2(mc.imageShifter(np.fft.fft2(G1), [-30, 10])))
        print(G1.shape, G2.shape)
        # plt.figure(0)
        # plt.imshow(np.concatenate((G1, G2)))
        # plt.show(block = True)
        imageCorr = mc.initial_correlation_image(np.fft.fft2(G1), np.fft.fft2(G2), 'phase', 2)

        # test to see if correlation image looks the same
        np.testing.assert_almost_equal(np.exp(1j * np.angle((np.multiply(np.fft.fft2(G2), np.conj(np.fft.fft2(G1)))))), imageCorr, decimal=4)

        # plt.figure(1)
        # plt.imshow(np.real(np.fft.ifft2(imageCorr)))
        # plt.show(block = True )
        out_phase = mc.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 2)
        self.assertEqual(out_phase, [-30.0, 10.0])

        # test cross correlation
        out_cross = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 1)
        self.assertEqual(out_cross, [-30.0, 10.0])

        # test hybrid correlation and dft upsample
        out_hybrid = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
        out_hybrid = list(out_hybrid)
        np.testing.assert_almost_equal(out_hybrid, [-30.0, 10.0], decimal = 6)

        # plt.imshow(np.subtract(G1, np.real(np.fft.ifft2(out.G2shift))))
        # plt.show(block = True)


    def test_different_shifts(self):
        '''
        Tests to check if mod is working correctly
        '''
        # filename = '/Users/Tom/Downloads/matt_beard.jpg'
        # G1 = cv2.imread(filename, 0)
        # G1 = G1[0:501, 0:501]
        G1 = np.zeros((101,101))
        G1[55,55] = 12
        # plt.imshow(G1)
        # plt.show(block = True)
        shifts = [[3., 1.], [-3., -1.], [-3., 1.], [3.,-1.]]
        shifts2 = [[10.3, 14.1], [-10.3, -14.1], [-10.3, 14.1], [10.3,-14.1]]
        for i in shifts:
            # print(i)
            G2 = np.real(np.fft.ifft2(mc.imageShifter(np.fft.fft2(G1), i)))
            # plt.imshow(np.concatenate((G1, G2)))
            # plt.show(block = True)
            with self.subTest(i = i):
                out_phase = mc.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 3)
                np.testing.assert_almost_equal(out_phase, i, decimal = 6)
                out_cross = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 3)
                np.testing.assert_almost_equal(out_cross, i, decimal = 6)
                out_hybrid = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
                out_hybrid = list(out_hybrid)
                np.testing.assert_almost_equal(out_hybrid, i, decimal = 6)
        for i in shifts2:
            # print(i)
            G2 = np.real(np.fft.ifft2(mc.imageShifter(np.fft.fft2(G1), i)))
            with self.subTest(i = i):

                out_phase = mc.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 10)
                out_phase = list(out_phase)
                np.testing.assert_almost_equal(out_phase, i, decimal = 2)

                out_hybrid = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 10)
                out_hybrid = list(out_hybrid)
                np.testing.assert_almost_equal(out_hybrid, i, decimal = 2)

                out_cross = mc.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 10)
                out_cross = list(out_cross)
                np.testing.assert_almost_equal(out_cross, i, decimal = 2)

if __name__ == "__main__":
    unittest.main()
