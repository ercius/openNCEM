"""
Tests for the basic functionality of the smv io module.
"""

import pytest

import time
from pathlib import Path
import tempfile
import numpy as np

import ncempy.io.smv


class Testsmv:
    """
    Test the SMV io module
    """

    @pytest.fixture
    def data_location(self):
        # Get the location of the test data files
        test_path = Path(__file__).resolve()
        root_path = test_path.parents[1]
        return root_path / Path('data')

    @pytest.fixture
    def temp_file(self):
        tt = tempfile.NamedTemporaryFile(mode='wb')
        tt.close()  # need to close the file to use it later
        return Path(tt.name)

    def test_write_read(self, temp_file):
        # Write out a temporary SMV file
        ncempy.io.smv.smvWriter(temp_file,
                                np.ones((10, 11), dtype=np.uint16))

        assert temp_file.exists() is True
        with ncempy.io.smv.fileSMV(temp_file) as f0:
            assert hasattr(f0, 'dataSize')
            assert f0.dataSize[0] == 10
    
    def test_biotin_file(self, data_location):
        file_path = data_location / Path('biotin_smv.img')
        with ncempy.io.smv.fileSMV(file_path) as f0:
            dd = f0.getDataset()
            assert 'data' in dd
            assert dd['data'].shape[0] == 2048