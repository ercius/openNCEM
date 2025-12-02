"""
Tests for the basic functionality of the dectris io module.
"""

import pytest

import time
from pathlib import Path
import tempfile
import numpy as np

import ncempy.io.dectris


class Testdectris:
    """
    Test the dectris io module
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    def test_read_data(self, data_location):
        file_path = data_location / Path('au_145mm_68kx_microprobe_01_master.h5')
        with ncempy.io.dectris.fileDECTRIS(file_path) as f0:
            dd = f0.getDataset()
            assert 'data' in dd

    def test_read_metadata(self, data_location):
        file_path = data_location / Path('au_145mm_68kx_microprobe_01_master.h5')
        with ncempy.io.dectris.fileDECTRIS(file_path) as f0:
            md = f0.getMetadata()
            if md:
                assert 'pixelSize' in md

    def test_str_input(self, data_location):
        file_path = data_location / Path('au_145mm_68kx_microprobe_01_master.h5')
        with ncempy.io.dectris.fileDECTRIS(str(file_path)) as f0:
            assert f0.data_shap[0] = 64
    
    #def test_read_file(self, data_location):
    #    file_path = data_location / Path('dectris_master.h5')
    #    with ncempy.io.dectris.fileDECTRIS(file_path) as f0:
    #        dd = f0.getDataset()
    #        assert 'data' in dd
    #        assert dd['data'].shape[0] == 2048
    
    # def test_dectrisReader(self, data_location):
    #     file_path = data_location / Path('dectris_master.h5')
    #     dd = ncempy.io.dectris.dectrisReader(file_path)
    #     assert 'data' in dd
    #     assert dd['data'].shape[0] == 2048
