from . import dm
from . import ser
from . import emd
from . import mrc
from . import emdVelox
from . import smv

from pathlib import Path


def read(filename, dsetNum=0):
    """
    Parameters
    ----------
    filename : str or pathlib.Path
        The path and name of the file to attempt to load. This chooses the Reader() function based on the file
        suffix.
    dsetNum : int
        The data set number to load if there are multiple data sets in a file. This is implemented for EMD and DM
        files only.
    """

    out = {}

    # check filename type
    if isinstance(filename, str):
        filename = Path(filename)
    elif isinstance(filename, Path):
        pass
    else:
        raise TypeError('Filename is supposed to be a string or pathlib.Path')

    if not filename.exists():
        raise FileNotFoundError

    # ensure lowercase suffix
    suffix = filename.suffix.lower()

    if suffix == '.ser':
        out = ser.serReader(filename)
    elif suffix in ('.dm3', '.dm4'):
        out = dm.dmReader(filename)
    elif suffix in ('.emd', '.h5', '.hdf5'):
        try:
            out = emd.emdReader(filename, dsetNum)
        except emd.NoEmdDataSets:
            out = emdVelox.emdVeloxReader(filename, dsetNum)
    elif suffix in ('.mrc', '.rec', '.st', '.ali'):
        out = mrc.mrcReader(filename)
    elif suffix in ('.smv', '.img'):
        out = smv.smvReader(filename)
    else:
        print('File suffix {} is not recognized.'.format(suffix))
        print('Supported formats are ser, mrc, rec, ali, st, emd, dm3, dm4, smv, img.')

    return out
