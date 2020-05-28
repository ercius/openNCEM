"""General tests for io module. More in depth testing is done in sepearte test_.py files for more advanced
capabilities.
"""

import pytest
import tempfile
from pathlib import Path

import numpy as np
import ncempy.io as nio

# Todo: Add emd velox (small file needed)

@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')

@pytest.fixture
def temp_file():
    tt = tempfile.NamedTemporaryFile(mode='wb')
    tt.close() # need to close the file to use it later
    return Path(tt.name)


def test_emd_berkeley(data_location):
    emd_path = data_location / Path('Acquisition_18.emd')
    with nio.emd.fileEMD(emd_path) as emd0:
        assert len(emd0.list_emds) == 1


def test_dm3_2d(data_location):
    dm_path = data_location / Path('dmTest_float32_nonSquare_diffPixelSize.dm3')
    with nio.dm.fileDM(dm_path) as dm3:
        dd = dm3.getDataset(0)
        assert dd['data'].shape == (64, 128)


def test_dm4_2d(data_location):
    dm_path = data_location / Path('dmTest_float32_nonSquare_diffPixelSize.dm4')
    with nio.dm.fileDM(dm_path) as dm4:
        dd = dm4.getDataset(0)
        assert dd['data'].shape == (64, 128)


def test_dm3_3d(data_location):
    dm_path = data_location / Path('dmTest_3D_int16_64,65,66.dm3')
    with nio.dm.fileDM(dm_path) as dm3:
        dd = dm3.getDataset(0)
        assert dd['data'].shape == (66, 65, 64)


def test_dm4_3d(data_location):
    dm_path = data_location / Path('dmTest_3D_int16_64,65,66.dm4')
    with nio.dm.fileDM(dm_path) as dm4:
        dd = dm4.getDataset(0)
        assert dd['data'].shape == (66, 65, 64)


def test_dm3_spectrum_1d(data_location):
    dm_path = data_location / Path('08_carbon.dm3')
    with nio.dm.fileDM(dm_path) as dm4:
        dd = dm4.getDataset(0)
        dd['data'].shape = 2048


def test_mrc(temp_file):

    # Write out a temporary mrc file
    nio.mrc.mrcWriter(temp_file,
                      np.ones((10, 11, 12), dtype=np.float32),
                      (1, 2, 3))

    assert temp_file.exists() is True

    with nio.mrc.fileMRC(temp_file) as mrc0:
        dd = mrc0.getDataset()
        assert dd['data'].shape == (10, 11, 12)