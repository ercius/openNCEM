import pytest

import numpy as np

from ncempy.algo import fourier_operations


def test_fourierShift():
    """Test on a simple dataset"""

    im = np.eye(10, 10)
    im_sh = fourier_operations.shiftImage(im, (2, 0))
    assert np.round(im_sh[2, 0]) == 1
    im_sh = fourier_operations.shiftImage(im, (0, 2))
    assert np.round(im_sh[0, 2]) == 1


def test_fourierRotate():
    """Test that it works"""
    im = np.eye(10, 10)

    theta = 3 * np.pi / 180.
    im_rot = fourier_operations.rotateImage(im, theta)

    im_rot = fourier_operations.rotateImage(im, theta, pad=True)
    assert im_rot.shape[0] == 12


def test_fourierShear():
    """Test that it works"""
    im = np.eye(10, 10)

    im_shear = fourier_operations.shearImage(im, 0, 0.5)

