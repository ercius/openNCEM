import pytest

from pathlib import Path

from scipy import ndimage
import numpy as np

import ncempy.eval as neval
import ncempy.io as nio

@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


def test_mutlicorr(data_location):

    with nio.emd.fileEMD(data_location / Path('Acquisition_18.emd')) as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    """
    Test to check if the correlation is working.
    """

    G1 = dd
    G2 = np.roll(dd, 10, axis=0)
    G2 = np.roll(G2, -5, axis=1)

    test = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G1), 'phase', 3)
    assert list(test) == [0, 0]

    out_phase = neval.multicorr(np.fft.fft2(G1), (np.fft.fft2(G2)), 'phase', 2)
    assert out_phase == [10.0, -5.0]

    # test cross correlation
    out_cross = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'cross', 1)
    assert out_cross == [10.0, -5.0]

    # test hybrid correlation and dft upsample
    out_hybrid = neval.multicorr(np.fft.fft2(G1), np.fft.fft2(G2), 'hybrid', 3)
    out_hybrid = list(out_hybrid)
    assert np.allclose(out_hybrid, [10, -5.0], rtol=1e-6)

def test_lineprofile():
    XX, YY = np.mgrid[0:100,0:100]
    #RR = np.sqrt((XX-50)**2 + (YY-50)**2)

    profile, (xx, yy) = neval.line_profile(XX, (0,0), (99, 99), 100)
    print(profile[0:5])
    assert np.allclose(profile[0:3],(0,1,2))

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
    out_stack, out_stack_shifts = neval.stack_align(stack0, align_type='static',
                                                    upsample_factor=up, method=method0)

    for sh_static, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_static[0] == -known_shift[0]
        assert sh_static[1] == -known_shift[1]

    # Test dynamic alignment (reference is n-1 image in stack)
    out_stack, out_stack_shifts = neval.stack_align(stack0, align_type='dynamic',
                                                    upsample_factor=up, method=method0)

    for sh_dynamic, known_shift in zip(shifts0, out_stack_shifts):
        assert sh_dynamic[0] == -known_shift[0]
        assert sh_dynamic[1] == -known_shift[1]
