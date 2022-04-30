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
