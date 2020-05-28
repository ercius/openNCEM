import pytest

from scipy import ndimage
import numpy as np

import ncempy.io as nio
from ncempy.eval.stack_align import stack_align


@pytest.fixture
def emd_path():
    return 'C:/Users/linol/Data/Acquisition_18.emd'


def test_stack_align(emd_path):
    up = 2
    method0 = 'hybrid'

    with nio.emd.fileEMD('C:/Users/linol/Data/Acquisition_18.emd') as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    stack0 = np.zeros((3, *dd.shape), dtype=dd.dtype)
    shifts0 = np.asarray(((range(0, -6, -2)), (range(0, 6, 2)))).T
    for ii, s in enumerate(shifts0):
        stack0[ii, :, :] = ndimage.shift(dd, s, mode='mirror')
    out_stack, out_stack_shifts = stack_align(stack0, align_type='static',
                                              upsample_factor=up, method=method0)
    #print('shifts0 = {}'.format(shifts0))
    #print('stack shifts = {}'.format(out_stack_shifts))
    for s0, s1 in zip(shifts0, out_stack_shifts):
        assert s0[0] == -s1[0]
        assert s0[1] == -s1[1]
