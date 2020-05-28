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

    # Todo: Update this to find the correct data directory
    @pytest.fixture
    def data_location(self):
        return Path(r'C:\Users\linol\scripting\openNCEMgh\ncempy\data')

    def test_read_ser_2d(self, data_location):

        # single 2D image file
        with ncempy.io.ser.fileSER(data_location / Path('16_STOimage_1.ser')) as ser0:
            assert ser0

            dd, md = ser0.getDataset(0)

            assert dd[0, 0] == 18024
            assert md['Calibration'][0]['CalibrationElement'] == 0

    def test_read_ser_3d(serl, data_location):
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
            emi = ser0.read_emi(str(data_location / Path('16_STOimage.emi')))

            assert emi['AcceleratingVoltage'] == 80000

    def test_write_emd(self):
        """
        Test the emd writing functionality.

        Todo: Update to work with pytest
        Todo: Use temporary directory and file to write out the converted EMD file.
        """

        # Make this test fail since its not updated yet
        assert True is False

        # single 2D image file
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser')
        if os.path.isfile('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd'):
            os.remove('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd')
        fser.writeEMD('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd')
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser', 'ncempy/test/resources/Pt_SAED_D910mm_single/im01.emi')
        if os.path.isfile('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd'):
            os.remove('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd')
        fser.writeEMD('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd')
        with self.assertRaises(IOError):
            fser.writeEMD('')
        
        # time series of 2D images
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01_1.ser', 'ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01.emi')
        if os.path.isfile('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd'):
            os.remove('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd')
        fser.writeEMD('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd')
