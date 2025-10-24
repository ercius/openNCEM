"""
Tests for the emdVelox io module.
"""

import pytest

from pathlib import Path

import ncempy.io.emdVelox


class TestEMDVelox:
    """
    Test the EMDVelox io module
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_readEMDVelox(self, data_location):
        dd0 = ncempy.io.emdVelox.emdVeloxReader(data_location / Path('STEM HAADF-DF4-DF2-BF Diffraction Micro.emd'),
                                                dsetNum=0)
        assert dd0['data'].ndim == 2

        dd2 = ncempy.io.emdVelox.emdVeloxReader(data_location / Path('STEM HAADF-DF4-DF2-BF Diffraction Micro.emd'),
                                                dsetNum=2)
        assert dd2['data'].ndim == 2

    def test_read_emd_stem(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('STEM HAADF Diffraction Micro.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2

        assert md['pixelSizeUnit'][0] == 'nm'
        assert md['pixelUnit'][0] == 'nm'
        assert round(md['pixelSize'][0], ndigits=4) == round(87.3708, ndigits=4)

        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('STEM HAADF-DF4-DF2-BF Diffraction Micro.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2

    def test_read_emd_tem(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('Camera Ceta Imaging Micro.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
            dd, md = emd0.get_dataset(emd0.list_emds[0])
        assert dd.ndim == 2

    def test_read_emd_diffraction(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('Camera Ceta Diffraction Micro.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2
        assert md['pixelSizeUnit'][0] == '1/m'
        assert md['pixelUnit'][0] == '1/m'

    def test_file_object(self, data_location):
        """Test fileEMDVelox class input with file object"""
        file_name = data_location / Path('STEM HAADF-DF4-DF2-BF Diffraction Micro.emd')
        fid = open(file_name, 'rb')
        emd0 = ncempy.io.emdVelox.fileEMDVelox(fid)
        assert hasattr(emd0, '_file_hdl')

    def test_metadata(self, data_location):
        file_path = data_location / Path('STEM HAADF Diffraction Micro.emd')
        with ncempy.io.emdVelox.fileEMDVelox(file_path) as emd0:
            md = emd0.getMetadata(0)
        assert md['AccelerationVoltage'] == '300000'
