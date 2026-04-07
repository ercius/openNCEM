"""
Tests for eels module.
"""

import pytest
import tempfile
from pathlib import Path

import numpy as np
import ncempy.io as nio
import ncempy.algo.eels as eels

@pytest.fixture
def data_location():
    # Get the location of the test data files
    test_path = Path(__file__).resolve()
    root_path = test_path.parents[1]
    return root_path / Path('data')


@pytest.fixture
def temp_file():
    tt = tempfile.NamedTemporaryFile(mode='wb')
    tt.close() # need to close the file to use it later
    return Path(tt.name)


def test_powerfit(data_location):
    data_path = data_location / Path('08_carbon.dm3')
    d0 = nio.read(data_path)

    bgnd, eLoss0 = eels.powerFit(d0['data'][0], [57, 392])

