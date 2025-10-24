# General tests for all of ncempy

import pytest
import tempfile
from pathlib import Path

import numpy as np
import ncempy


@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


@pytest.fixture
def temp_file():
    tt = tempfile.NamedTemporaryFile(mode='wb')
    tt.close()  # need to close the file to use it later
    return Path(tt.name)


def test_read(data_location):
    emd_path = data_location / Path('Acquisition_18.emd')
    _ = ncempy.read(emd_path)
