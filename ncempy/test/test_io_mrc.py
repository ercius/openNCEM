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

    @staticmethod
    def _add_extended_header(path, ext_header, ext_type=b'\x00' * 4):
        # Insert an extended header after the 1024 byte main header
        raw = bytearray(path.read_bytes())
        raw[92:96] = np.int32(len(ext_header)).tobytes()  # NSYMBT
        raw[104:108] = ext_type  # EXTTYP
        path.write_bytes(bytes(raw[:1024]) + ext_header + bytes(raw[1024:]))

    def test_symmetry_extended_header(self, temp_file):
        # Two 80 character symmetry records, as in older EMDB maps without EXTTYP
        ncempy.io.mrc.mrcWriter(temp_file, np.ones((10, 11, 12), dtype=np.float32), (1, 2, 3))
        self._add_extended_header(temp_file, b'X,  Y,  Z'.ljust(80) * 2)

        with ncempy.io.mrc.fileMRC(temp_file) as f0:
            assert f0.FEIinfo == {}
            assert np.allclose(f0.voxelSize, (1, 2, 3))
            assert f0.getDataset()['data'].shape == (10, 11, 12)

    def test_fei1_extended_header(self, temp_file):
        # FEI1 metadata blocks are not the legacy FEI layout
        ncempy.io.mrc.mrcWriter(temp_file, np.ones((10, 11, 12), dtype=np.float32), (1, 2, 3))
        self._add_extended_header(temp_file, np.full(10 * 768 // 4, 1e-5, dtype=np.float32).tobytes(), b'FEI1')

        with ncempy.io.mrc.fileMRC(temp_file) as f0:
            assert f0.FEIinfo == {}
            assert np.allclose(f0.voxelSize, (1, 2, 3))

    def test_legacy_fei_extended_header(self, temp_file):
        # Legacy FEI: 32 floats per section, pixel size in meters is the 12th
        ncempy.io.mrc.mrcWriter(temp_file, np.ones((10, 11, 12), dtype=np.float32), (1, 2, 3))
        ext = np.zeros((1024, 32), dtype=np.float32)
        ext[:, 11] = 2.5e-10
        self._add_extended_header(temp_file, ext.tobytes())

        with ncempy.io.mrc.fileMRC(temp_file) as f0:
            assert f0.FEIinfo['pixel_size'] == pytest.approx(2.5e-10)
            assert np.allclose(f0.voxelSize, (1, 2.5, 2.5))

    # def test_metadata(self):
    #     file_path = Path('/mnt/NAS-NCEM_Data/TitanX/KateG/KateG/Greenlee_20190413_FeFeO_KGBox4G2/tiltSeries_20190413_FeFeO_neg70to65.mrc')
    #     with ncempy.io.mrc.fileMRC(file_path) as f0:
    #         md = f0.getMetadata()
    #     md['tilt_axis'])
