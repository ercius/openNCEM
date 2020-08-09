"""
Tests for the basic functionality of the dm io module.
"""

import pytest

import time
import inspect
import os
from pathlib import Path

import ncempy.io.dm


class Testdm3:
    """
    Test the DM3 io module
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def _read_dm3_data(self, file_path, on_memory=False):
        """Creates a fileDM and reads its data metadata

        Parameters
        ----------
            file_path : pathlib.path
                The full path of the file to load

            on_memory : bool
                if True, the dm file will be opened in on memory mode.
        """

        with ncempy.io.dm.fileDM(file_path, on_memory=on_memory) as f:
            if on_memory:
                assert f._on_memory
            f.parseHeader()
            ds = f.getDataset(0)
            img = ds['data']

            if img.ndim > 4:
                raise ValueError("Images with more than four dimensions not supported.")
            metadata = dict(dimensions=img.ndim, header=f.allTags,
                            metadata={x: ds[x] for x in ds.keys() if x != "data"})

        return metadata, img

    def test_read_dm3(self, data_location):
        
        metadata, img = self._read_dm3_data(data_location / Path('dmTest_3D_int16_64,65,66.dm3'))

        assert metadata["dimensions"] == 3
        assert metadata["metadata"]["pixelOrigin"] == [0.0, 0.0, 0.0]
        assert img.ndim == 3

    def test_read_dm3_1d(self, data_location):
        
        metadata, img = self._read_dm3_data(data_location / Path('08_carbon.dm3'))

        assert metadata["dimensions"] == 2
        assert metadata["metadata"]["pixelOrigin"] == [0.0, -2400.0]
        assert img.shape == (1, 2048)

    def test_read_dm3_on_memory(self, data_location):

        metadata, img = self._read_dm3_data(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=True)

        assert metadata["dimensions"] == 3
        assert metadata["metadata"]["pixelOrigin"] == [0.0, 0.0, 0.0]

    def test_dm4_memory_vs_file_performance(self, data_location):
        """ Test speed improvement with on_memory option.
        Even with a local HDD, memory read should be x10 faster.

        """
        m0 = time.time()
        for i in range(10):
            _ = self._read_dm3_data(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=False)
            delta0 = time.time() - m0
            
            m1 = time.time()
            _ = self._read_dm3_data(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=True)

            delta1 = time.time() - m1

            assert delta0 > delta1

    def test_extract_on_memory(self, data_location):
        fPath = data_location / Path('dmTest_3D_int16_64,65,66.dm3')
        with ncempy.io.dm.fileDM(fPath, on_memory=False) as f:
            ds = f.getDataset(0)
            img3d_no_on_mem = ds['data']

        with ncempy.io.dm.fileDM(fPath, on_memory=True) as f:
            ds = f.getDataset(0)
            img3d_on_mem = ds['data']

        assert img3d_on_mem[0, 0, 0] == img3d_no_on_mem[0, 0, 0]

    def test_on_memory(self, data_location):
        """ Ensure that data loaded is available in memory after the file is closed. This test might crash
        python if the memory is no longer available. See also test_dmReader
        """
        with ncempy.io.dm.fileDM(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=True) as f:
            dv = f.getDataset(0)
            ds = f.getSlice(0, 0)

            # Test Data is in memory with file open
            assert dv['data'][0, 0, 0] == 0
            assert ds['data'][0, 0] == 0

        # Ensure it is still available after closing the file
        assert dv['data'][0, 0, 0] == 0
        assert ds['data'][0, 0] == 0

    def test_memmap(self, data_location):
        import numpy as np
        with ncempy.io.dm.fileDM(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=True) as f:
            m = f.getMemmap(0)
            assert isinstance(m, np.ndarray)
            assert isinstance(m, np.memmap)
            assert m[0, 0, 0] == 0  # test while file open

        assert m[0, 0, 0] == 0  # test when file closed
        del m  # close the memmap

        with ncempy.io.dm.fileDM(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=False) as f:
            m = f.getMemmap(0)
            assert isinstance(m, np.ndarray)
            assert isinstance(m, np.memmap)
            assert m[0, 0, 0] == 0

        assert m[0, 0, 0] == 0
        del m

        with ncempy.io.dm.fileDM(data_location / Path('08_carbon.dm3'), on_memory=False) as f:
            m = f.getMemmap(0)
            assert isinstance(m, np.ndarray)
            assert isinstance(m, np.memmap)
            assert int(m[0]) == 21281

    def test_dmReader(self, data_location):
        """Test that the simplified dmReader function works and loads the data into memory. If test_on_memory
        fails then this will likely fail for the same reason."""
        dm0 = ncempy.io.dm.dmReader(data_location / Path('dmTest_3D_int16_64,65,66.dm3'))

        assert dm0['data'][0, 0, 0] == 0

        dm1 = ncempy.io.dm.dmReader(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=True)
        assert dm1['data'][0, 0, 0] == 0
        dm2 = ncempy.io.dm.dmReader(data_location / Path('dmTest_3D_int16_64,65,66.dm3'), on_memory=False)
        assert dm2['data'][0, 0, 0] == 0

    def test_writeTags(self, data_location):
        file_name = data_location / Path('08_carbon.dm3')
        with ncempy.io.dm.fileDM(file_name) as dm0:
            dm0.writeTags()
            new_loc = data_location / Path(file_name.stem + '_tags.txt')
            assert new_loc.exists()
            new_loc.unlink()
            dm0.writeTags(new_folder_path_for_tags=data_location)
            assert new_loc.exists()
            new_loc.unlink()

    def test_dmReader_spectra(self, data_location):
        # Ensure that the coordinates in DM spectra files are correct
        file_name = data_location / Path('08_carbon.dm3')
        dm0 = ncempy.io.dm.dmReader(file_name)
        assert dm0['data'].ndim == 2
        assert round(dm0['coords'][1][0]) == 240
        assert round(dm0['coords'][1][-1], ndigits=1) == round(444.7, ndigits=1)
