"""
Tests for the eval.multicorr module and eval.multicorr_funcs
"""

import pytest

import numpy as np
import ncempy.eval as neval
import ncempy.algo as nalgo


class TestShiftedimages():
    """Test the multicorr module on two shifted images"""

    def test_multicorr(self):
        """
        Test to check if the correlation is working.
        """
        g2 = np.zeros((3, 3), dtype=np.float)
        test = neval.multicorr(np.fft.fft2(g2), np.fft.fft2(g2), 'phase', 3)
        assert list(test) == [0, 0]

        G1 = np.zeros((100, 100))
        G1[42, 35] = 12
        G2 = np.real(np.fft.ifft2(nalgo.multicorr_funcs.imageShifter(np.fft.fft2(G1), [-30, 10])))
        imageCorr = nalgo.multicorr_funcs.initial_correlation_image(np.fft.fft2(G1), np.fft.fft2(G2), 'phase', 2)
        
        # test to see if correlation image looks the same
        assert np.allclose(np.exp(1j * np.angle((np.multiply(np.fft.fft2(G2), np.conj(np.fft.fft2(G1)))))), imageCorr,
                           rtol=1e-6)

        out_phase = neval.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 2)
        assert out_phase == [-30.0, 10.0]

        # test cross correlation
        out_cross = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 1)
        assert out_cross == [-30.0, 10.0]

        # test hybrid correlation and dft upsample
        out_hybrid = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
        out_hybrid = list(out_hybrid)
        assert np.allclose(out_hybrid, [-30, 10.0], rtol=1e-6)

    def test_different_shifts(self):
        """
        Tests to check if mod is working correctly
        """
        G1 = np.zeros((101, 101))
        G1[55, 55] = 12
        shifts = [[3., 1.], [-3., -1.], [-3., 1.], [3.,-1.]]
        shifts2 = [[10.3, 14.1], [-10.3, -14.1], [-10.3, 14.1], [10.3, -14.1]]
        for i in shifts:
            G2 = np.real(np.fft.ifft2(nalgo.multicorr_funcs.imageShifter(np.fft.fft2(G1), i)))
            out_phase = neval.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 3)
            assert np.allclose(out_phase, i, rtol=1e-6)
            out_cross = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 3)
            assert np.allclose(out_phase, i, rtol=1e-6)
            out_hybrid = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
            out_hybrid = list(out_hybrid)
            assert np.allclose(out_phase, i, rtol=1e-6)

        for i in shifts2:
            G2 = np.real(np.fft.ifft2(nalgo.multicorr_funcs.imageShifter(np.fft.fft2(G1), i)))

            out_phase = neval.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 10)
            out_phase = list(out_phase)
            assert np.allclose(out_phase, i, rtol=1e-6)

            out_hybrid = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 10)
            out_hybrid = list(out_hybrid)
            assert np.allclose(out_phase, i, rtol=1e-6)

            out_cross = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 10)
            out_cross = list(out_cross)
            assert np.allclose(out_phase, i, rtol=1e-6)
