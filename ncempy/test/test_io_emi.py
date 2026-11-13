"""
Tests for the emi io module.
"""

import pytest

from pathlib import Path

import ncempy.io.emi


class TestEmi():
    """
    Test the EMI io module.
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_read_emi_image(self, data_location):
        with ncempy.io.emi.fileEMI(data_location / Path('flucam_image.emi')) as emi0:
            assert emi0

            emi0.parse_file()

            assert len(emi0.image_locations) == 1
            assert emi0.image_name == ['FluCam CCD']
            assert emi0.image_size == [(512, 512)]
            assert emi0.data_type == ['<u4']

            dd = emi0.getDataset(0)

            assert dd['name'] == 'FluCam CCD'
            assert dd['data'].shape == (512, 512)
            assert dd['data'].dtype == '<u4'
            assert dd['data'][256, 256] == 236
            assert dd['data'].sum() == 3044561

    def test_getDataset_bad_index(self, data_location):
        with ncempy.io.emi.fileEMI(data_location / Path('flucam_image.emi')) as emi0:
            emi0.parse_file()

            with pytest.raises(IndexError):
                emi0.getDataset(999)
