import pytest

from ncempy.algo import rebin

import numpy as np


def test_rebin():
    aa = np.ones((50, 50), dtype='<u2')

    # Test sum and uint16
    bb = rebin(aa, 2, funcType='sum')

    assert bb[0, 0] == 4

    aa = np.ones((50, 50), dtype='f')

    # Test mean and float
    bb = rebin(aa, 2, funcType='mean')

    assert bb[0, 0] == 1.0
