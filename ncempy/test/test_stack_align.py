import pytest

from pathlib import Path

from scipy import ndimage
import numpy as np

import ncempy.io as nio
from ncempy.eval.stack_align import stack_align


@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


def test_stack_align(data_location):
    up = 2
    method0 = 'hybrid'

    with nio.emd.fileEMD(data_location / Path('Acquisition_18.emd')) as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    stack0 = np.zeros((3, *dd.shape), dtype=dd.dtype)
    shifts0 = np.asarray(((range(0, -6, -2)), (range(0, 6, 2)))).T
    for ii, s in enumerate(shifts0):
        stack0[ii, :, :] = ndimage.shift(dd, s, mode='mirror')
    out_stack, out_stack_shifts = stack_align(stack0, align_type='static',
                                              upsample_factor=up, method=method0)

    for s0, s1 in zip(shifts0, out_stack_shifts):
        assert s0[0] == -s1[0]
        assert s0[1] == -s1[1]
