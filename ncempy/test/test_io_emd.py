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
        tt.close()  # need to close the file to use it later
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
        
        try:
            with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:
                emd0.put_comment('file created')
        except ncempy.io.emd.NoEmdDataSets:
            pass

        try:
            with ncempy.io.emd.fileEMD(temp_file, readonly=True) as emd1:
                assert len(emd1.comments.attrs) == 1
        except ncempy.io.emd.NoEmdDataSets:
            pass

    def test_write_emd(self, temp_file):
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        dims = ncempy.io.emd.defaultDims(dd, pixel_size=(0.1, 0.2,  0.3), pixel_unit=('nm','nm','nm'))

        with ncempy.io.emd.fileEMD(temp_file, readonly=False) as emd0:
            emd0.put_emdgroup('pytest', dd, dims)

        with ncempy.io.emd.fileEMD(temp_file, readonly=True) as emd1:
            dd2, dims2 = emd1.get_emdgroup(0)
        assert dd2.shape == (10, 11, 12)
        assert np.allclose(dims2[0][0][1], 0.1)
        assert np.allclose(dims2[1][0][1], 0.2)
        assert np.allclose(dims2[2][0][1], 0.3)
        assert dims2[0][2] == 'nm'
        assert dims2[1][2] == 'nm'
        assert dims2[2][2] == 'nm'
    
    def test_emd_writer_simple(self, temp_file):
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        ncempy.io.emd.emdWriter(temp_file, dd)
    
    def test_emd_writer(self, temp_file):
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        ncempy.io.emd.emdWriter(temp_file, dd, pixel_size=(1, 2, 3), pixel_unit=('nm','Ang','um'))
        with ncempy.io.emd.fileEMD(temp_file) as emd0:
            dd2, dims2 = emd0.get_emdgroup(0)
            assert dd2.ndim == 3
            assert dims2[0][0].size == 10
            assert dims2[1][0].size == 11
            assert dims2[2][0].size == 12
            # Test pixel size is stored correctly
            assert int(dims2[0][0][1] - dims2[0][0][0]) == 1
            assert int(dims2[1][0][1] - dims2[1][0][0]) == 2
            assert int(dims2[2][0][1] - dims2[2][0][0]) == 3
            assert dims2[0][2] == 'nm'
            assert dims2[1][2] == 'Ang'
            assert dims2[2][2] == 'um'

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

    def test_file_object(self, data_location):
        # Test fileEMD class input with file object
        file_name = data_location / Path('Acquisition_18.emd')
        fid = open(file_name, 'rb')
        emd0 = ncempy.io.emd.fileEMD(fid)
        assert hasattr(emd0, 'file_hdl')

    def test_memmap(self, data_location):
        emd1 = ncempy.io.emd.fileEMD(data_location / Path('Acquisition_18.emd'))
        d, dims = emd1.get_memmap(0)
        del emd1
        assert d[0, 0] == 12487

    def test_memmap_in_function(self, data_location):
        f = data_location / Path('Acquisition_18.emd')

        def load_memmap(fpath, n):
            f0 = ncempy.io.emd.fileEMD(fpath)
            data0, dims0 = f0.get_memmap(n)
            return data0, dims0

        data, dims = load_memmap(f, 0)
        assert data[0, 0] == 12487

    def test_bad_dims(self, temp_file):
        import h5py
        # Create a data set with missing attributes in the dim vectors
        with h5py.File(temp_file, 'w') as f0:
            gg = f0.create_group('/data/temp')
            gg.attrs['emd_group_type'] = 1
            gg.create_dataset('data', data=np.zeros((10, 10)))
            dim1 = gg.create_dataset('dim1', data=(0, 1))
            # dim1.attrs['name'] = 'X'
            dim1.attrs['units'] = 'n_m'
            dim2 = gg.create_dataset('dim2', data=(0, 1))
            dim2.attrs['name'] = 'Y'
            # dim2.attrs['units'] = 'n_m'

        with ncempy.io.emd.fileEMD(temp_file) as f0:
            _, dims = f0.get_emdgroup(0)
        assert dims[0][1] == 'dim1'
        assert dims[1][2] == 'pixels'

    def test_version(self, temp_file):
        import h5py
        dd = np.ones((10, 11, 12), dtype=np.uint16)
        ncempy.io.emd.emdWriter(temp_file, dd)
        # Change the minor version manually
        with h5py.File(temp_file, 'r+') as f1:
            f1.attrs['version_minor'] = 3
        try:
            dd = ncempy.io.emd.emdReader(temp_file)
        except ncempy.io.emd.UnsupportedEmdVersion:
            pass
        
        # This should not error since an EMD daat set was found
        dd = ncempy.io.emd.emdReader(temp_file)
        
    def test_no_emds(self, temp_file, data_location):
        import h5py
        # Creat a HDF5 file with a data set
        with h5py.File(temp_file, 'w') as f0:
            f0.create_dataset('test', data=np.zeros((10,10)))
            
        # Try to read the data with emdReader
        try:
            ncempy.io.emd.emdReader(temp_file)
        except ncempy.io.emd.NoEmdDataSets:
            pass
        
        # Try to read a Velox file
        try:
            ncempy.io.emd.emdReader(data_location / Path('STEM HAADF-DF4-DF2-BF Diffraction Micro.emd'))
        except ncempy.io.emd.NoEmdDataSets:
            pass
