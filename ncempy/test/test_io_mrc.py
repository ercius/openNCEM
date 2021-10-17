"""
Tests for the basic functionality of the mrc io module.
"""

import pytest

from pathlib import Path
import tempfile
import numpy as np

import ncempy.io.mrc


class Testmrc:
    """
    Test the DM3 io module
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

    def test_file_object(self, temp_file):
        # Write out a temporary mrc file
        ncempy.io.mrc.mrcWriter(temp_file,
                                np.ones((10, 11, 12), dtype=np.float32),
                                (1, 2, 3))

        assert temp_file.exists() is True
        fid = open(temp_file, 'rb')
        dm0 = ncempy.io.mrc.fileMRC(fid)
        assert hasattr(dm0, 'fid')
