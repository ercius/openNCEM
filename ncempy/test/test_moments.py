import pytest

import numpy as np
from scipy import ndimage

from ncempy.algo import moments


def test_moments():
    """Test on a rotated line"""

    ang = 35  # degrees
    sh = (2, -2)

    # Create a a simple line
    aa = np.zeros((101, 101))
    aa[50, 1:-1] = 1
    aa = ndimage.rotate(aa, ang, order=3, reshape=False)
    aa = ndimage.shift(aa, sh)

    M = moments.moments(aa)
    mc = moments.moments_central(aa)

    c = moments.centroid(M)  # the center
    print('centroid = {}'.format(c))
    for mm, nn in zip(sh, c):
        assert mm + 50 == round(nn)

    th = moments.moment_angle(mc) + np.pi/2.  # in radians
    assert ang == round(th * 180 / np.pi)
