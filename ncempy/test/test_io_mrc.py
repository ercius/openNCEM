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
    Test the MRC io module
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
        with open(temp_file, 'rb') as f0:
            mrc0 = ncempy.io.mrc.fileMRC(f0)
            assert hasattr(mrc0, 'fid')

    def test_min_max_mean(self, temp_file):
        # DMIN, DMAX and DMEAN are floats in the MRC header
        data = np.arange(10 * 11 * 12, dtype=np.float32).reshape((10, 11, 12)) - 100.5
        ncempy.io.mrc.mrcWriter(temp_file, data, (1, 2, 3))

        with ncempy.io.mrc.fileMRC(temp_file) as f0:
            assert f0.minMaxMean.dtype == np.float32
            assert f0.minMaxMean[0] == pytest.approx(data.min())
            assert f0.minMaxMean[1] == pytest.approx(data.max())
            assert f0.minMaxMean[2] == pytest.approx(data.mean())

    # def test_metadata(self):
    #     file_path = Path('/mnt/NAS-NCEM_Data/TitanX/KateG/KateG/Greenlee_20190413_FeFeO_KGBox4G2/tiltSeries_20190413_FeFeO_neg70to65.mrc')
    #     with ncempy.io.mrc.fileMRC(file_path) as f0:
    #         md = f0.getMetadata()
    #     md['tilt_axis'])
