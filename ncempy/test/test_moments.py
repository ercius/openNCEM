import pytest

import numpy as np

from ncempy.algo import moments
from ncempy.algo import fourier_operations


def test_moments():
    """Test on a rotated line"""

    ang = 35  # degree
    sh = (2, -2)

    # Create a a simple line
    aa = np.zeros((101, 101))
    aa[50, 1:-1] = 1
    aa = fourier_operations.rotateImage(aa, ang * np.pi / 180.)
    aa = fourier_operations.shiftImage(aa, sh)

    M = moments.moments(aa)
    mc = moments.moments_central(aa)

    c = moments.centroid(M)  # the center
    print('centroid = {}'.format(c))
    for mm, nn in zip(sh, c):
        assert mm + 50 == round(nn)

    th = moments.moment_angle(mc) + np.pi/2.  # in radians
    assert ang == round(th * 180 / np.pi)


def test_moments_dtype():
    """Test that moments are calculated correctly with uint8 input dtype"""
    # Create test image
    YY, XX = np.meshgrid(np.linspace(0, 24, 25), np.linspace(0, 24, 25))
    RR = np.sqrt((XX - 24)**2 + (YY - 24)**2)
    RR = 255 * np.abs((RR / RR.max()) - 1)

    c1 = moments.centroid(moments.moments(RR))
    c2 = moments.centroid(moments.moments(RR.astype(np.uint8)))

    assert int(c1[0]) == int(c2[0])
    assert int(c1[1]) == int(c2[1])
