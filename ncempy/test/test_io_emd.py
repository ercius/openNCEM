"""
Tests for the emd io module.
"""

import pytest

from pathlib import Path
import tempfile

import numpy as np

import ncempy.io.emd


class Testemd:
    """
    Test the EMD io module
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    @pytest.fixture
    def temp_location(self):
        tmp_dir = tempfile.TemporaryDirectory()
        return Path(tmp_dir.name)

    @pytest.fixture
    def temp_file(self):
        tt = tempfile.NamedTemporaryFile(mode='wb')
        tt.close() # need to close the file to use it later
        return Path(tt.name)

    def test_read_emd_2d(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emd.fileEMD(data_location / Path('Acquisition_18.emd'), readonly=True) as emd0:
            dd, md = emd0.get_emdgroup(emd0.list_emds[0])
            assert dd.ndim == 2

            dd, md = emd0.get_emdgroup(0)
            assert dd.ndim == 2
        
    def test_comments(self, temp_file):
        """Test putting and retrieving comments"""

        with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:
            emd0.put_comment('file created')

        with ncempy.io.emd.fileEMD(temp_file, readonly=True) as emd1:
            assert len(emd1.comments.attrs) == 1

    def test_write_emd(self, temp_file):
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        dims = ncempy.io.emd.defaultDims(dd, pixel_size=(0.1, 0.2,  0.3))

        with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:
            emd0.put_emdgroup('pytest', dd, dims)

        with ncempy.io.emd.fileEMD(temp_file, readonly=True) as emd1:
            dd2, dims2 = emd1.get_emdgroup(0)
        assert dd2.shape == (10, 11, 12)

    def test_emd_writer(self, temp_file):
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        ncempy.io.emd.emdWriter(temp_file, dd, pixel_size=(1, 2, 3))
        with ncempy.io.emd.fileEMD(temp_file) as emd0:
            dd2, dims2 = emd0.get_emdgroup(0)
            assert dd2.ndim == 3
            assert dims2[0][0].size == 10
            assert dims2[1][0].size == 11
            assert dims2[2][0].size == 12
            # Test pixel size is stored correctly
            assert int(dims2[2][0][1] - dims2[2][0][0]) == 3

    def test_dim_string(self, data_location):
        """Test for dims with string attributes as name and units"""
        with ncempy.io.emd.fileEMD(data_location / Path('emd_type1_stringDims.h5')) as emd0:
            dims = emd0.get_emddims(emd0.list_emds[0])
        assert len(dims) == 4

        out = ncempy.io.emd.emdReader(data_location / Path('emd_type1_stringDims.h5'))
        assert out['data'].ndim == 4

    def test_file_exists(self, temp_file):
        dd = np.zeros((10, 11, 12))
        ncempy.io.emd.emdWriter(temp_file, dd, pixel_size=(1, 2, 3))
        assert temp_file.exists()
        try:
            ncempy.io.emd.emdWriter(temp_file, dd, pixel_size=(1, 2, 3), overwrite=True)
        except FileExistsError:
            assert False
