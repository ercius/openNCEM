import pytest

from pathlib import Path

from scipy import ndimage
import numpy as np

import ncempy.io as nio
from ncempy.eval.stack_align import stack_align
from ncempy.algo.align import stack_align as stack_align_algo


@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


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
