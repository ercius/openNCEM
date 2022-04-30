import pytest

from pathlib import Path

from scipy import ndimage
import numpy as np

import ncempy.algo as nalgo
import ncempy.io as nio
from ncempy.eval.stack_align import stack_align
from ncempy.algo import stack_align as stack_align_algo


@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


def test_rebin():
    aa = np.ones((50, 50), dtype='<u2')

    # Test sum and uint16
    bb = nalgo.rebin(aa, 2, funcType='sum')

    assert bb[0, 0] == 4

    aa = np.ones((50, 50), dtype='f')

    # Test mean and float
    bb = nalgo.rebin(aa, 2, funcType='mean')

    assert bb[0, 0] == 1.0

def test_stack_align_multicorr(data_location):
    up = 2
    method0 = 'hybrid'

    with nio.emd.fileEMD(data_location / Path('Acquisition_18.emd')) as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    stack0 = np.zeros((3, *dd.shape), dtype=dd.dtype)
    shifts0 = np.asarray(((range(0, -6, -2)), (range(0, 6, 2)))).T
    for ii, s in enumerate(shifts0):
        stack0[ii, :, :] = ndimage.shift(dd, s, mode='mirror')

    # Test static alignment (reference is first image in stack)
    out_stack, out_stack_shifts = stack_align(stack0, align_type='static',
                                              upsample_factor=up, method=method0)

    for sh_static, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_static[0] == -known_shift[0]
        assert sh_static[1] == -known_shift[1]

    # Test dynamic alignment (reference is n-1 image in stack)
    out_stack, out_stack_shifts = stack_align(stack0, align_type='dynamic',
                                              upsample_factor=up, method=method0)

    for sh_dynamic, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_dynamic[0] == -known_shift[0]
        assert sh_dynamic[1] == -known_shift[1]


def test_stack_align_cross(data_location):

    with nio.emd.fileEMD(data_location / Path('Acquisition_18.emd')) as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    # Create a stack with known shifts
    stack0 = np.zeros((3, *dd.shape), dtype=dd.dtype)
    shifts0 = np.asarray(((range(0, -6, -2)), (range(0, 6, 2)))).T
    for ii, s in enumerate(shifts0):
        stack0[ii, :, :] = ndimage.shift(dd, s, mode='mirror')

    # Test static alignment (reference is first image in stack)
    out_stack, out_stack_shifts = stack_align_algo(stack0, align_type='static')

    for sh_static, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_static[0] == -known_shift[0]
        assert sh_static[1] == -known_shift[1]

    # Test dynamic alignment (reference is n-1 image in stack)
    out_stack, out_stack_shifts = stack_align_algo(stack0, align_type='dynamic')

    for sh_dynamic, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_dynamic[0] == -known_shift[0]
        assert sh_dynamic[1] == -known_shift[1]


def test_moments():
    """Test on a rotated line"""

    ang = 35  # degree
    sh = (2, -2)

    # Create a a simple line
    aa = np.zeros((101, 101))
    aa[50, 1:-1] = 1
    aa = nalgo.rotateImage(aa, ang * np.pi / 180.)
    aa = nalgo.shiftImage(aa, sh)

    M = nalgo.moments(aa)
    mc = nalgo.moments_central(aa)

    c = nalgo.centroid(M)  # the center
    print('centroid = {}'.format(c))
    for mm, nn in zip(sh, c):
        assert mm + 50 == round(nn)

    th = nalgo.moment_angle(mc) + np.pi/2.  # in radians
    assert ang == round(th * 180 / np.pi)


def test_moments_dtype():
    """Test that moments are calculated correctly with uint8 input dtype"""
    # Create test image
    YY, XX = np.meshgrid(np.linspace(0, 24, 25), np.linspace(0, 24, 25))
    RR = np.sqrt((XX - 24)**2 + (YY - 24)**2)
    RR = 255 * np.abs((RR / RR.max()) - 1)

    c1 = nalgo.centroid(nalgo.moments(RR))
    c2 = nalgo.centroid(nalgo.moments(RR.astype(np.uint8)))

    assert int(c1[0]) == int(c2[0])
    assert int(c1[1]) == int(c2[1])


def test_fourierShift():
    """Test on a simple dataset"""

    im = np.eye(10, 10)
    im_sh = nalgo.shiftImage(im, (2, 0))
    assert np.round(im_sh[2, 0]) == 1
    im_sh = nalgo.shiftImage(im, (0, 2))
    assert np.round(im_sh[0, 2]) == 1


def test_fourierRotate():
    """Test that it works"""
    im = np.eye(10, 10)

    theta = 3 * np.pi / 180.
    im_rot = nalgo.rotateImage(im, theta)

    im_rot = nalgo.rotateImage(im, theta, pad=True)
    assert im_rot.shape[0] == 12


def test_fourierShear():
    """Test that it works"""
    im = np.eye(10, 10)

    im_shear = nalgo.shearImage(im, 0, 0.5)
