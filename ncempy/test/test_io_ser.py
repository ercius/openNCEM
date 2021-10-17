"""
Tests for the ser io module.
"""

import pytest

from pathlib import Path
import os
import os.path

import numpy as np

import ncempy.io.ser


class Testser():
    """
    Test the SER io module.
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_read_ser_2d(self, data_location):

        # single 2D image file
        with ncempy.io.ser.fileSER(data_location / Path('16_STOimage_1.ser')) as ser0:
            assert ser0

            dd, md = ser0.getDataset(0)

            assert dd[0, 0] == 18024
            assert md['Calibration'][0]['CalibrationElement'] == 0

    def test_read_ser_3d(self, data_location):
        with ncempy.io.ser.fileSER(data_location / Path('01_Si110_5images_1.ser')) as ser0:
            assert ser0

            dd, md = ser0.getDataset(1)

            assert dd[0, 0] == 14167

    def test_read_emi(self, data_location):
        """
        Test the emi reading functionality.
        """
        
        # ser file for this testing
        with ncempy.io.ser.fileSER(data_location / Path('16_STOimage_1.ser')) as ser0:
            assert ser0._emi['AcceleratingVoltage'] == 80000

        # Test reading emi file without ser file
        emi = ncempy.io.ser.read_emi(str(data_location / Path('16_STOimage.emi')))
        assert emi['AcceleratingVoltage'] == 80000

    def test_serReader(self, data_location):
        dd = ncempy.io.ser.serReader(data_location / Path('16_STOimage_1.ser'))

        assert dd['data'][0, 0] == 18024

    def test_file_object(self, data_location):
        # Test fileSER class input with file object
        file_name = data_location / Path('16_STOimage_1.ser')
        fid = open(file_name, 'rb')
        ser0 = ncempy.io.ser.fileSER(fid)
        assert hasattr(ser0, '_file_hdl')
