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
        #test_path = Path(__file__).resolve()
        #root_path = test_path.parents[1]
        root_path = Path(r'C:\Users\linol')
        return root_path / Path('data')

    def test_read_emd_stem(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('1435 1.2 Mx STEM HAADF-DF4-DF2-BF.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2

        assert md['pixelSizeUnit'][0] == 'nm'
        assert round(md['pixelSize'][0], ndigits=4) == round(0.0849,ndigits=4)

    def test_read_emd_tem(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('1532 Ceta 382.3 kx Imaging.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2

    def test_read_emd_diffraction(self, data_location):
        """Test loading data with list_emds and just an integer"""
        with ncempy.io.emdVelox.fileEMDVelox(data_location / Path('1534 Ceta 378 mm Diffraction.emd')) as emd0:
            dd, md = emd0.get_dataset(0)
        assert dd.ndim == 2
        assert md['pixelSizeUnit'][0] == '1/m'