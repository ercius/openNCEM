import ncempy.io as nio

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
        out = nio.ser.serReader(filename)
    elif suffix in ('.dm3', '.dm4'):
        out = nio.dm.dmReader(filename)
    elif suffix in ('.emd', '.h5', '.hdf5'):
        is_velox = False
        with nio.emd.fileEMD(filename) as emd0:
            if len(emd0.list_emds) > 0:
                out = nio.emd.emdReader(filename, dsetNum)
            else:
                is_velox = True
        if is_velox:
            out = nio.emdVelox.emdVeloxReader(filename, dsetNum)
    elif suffix in ('.mrc', '.rec', '.st', '.ali'):
        out = nio.mrc.mrcReader(filename)
    else:
        print('File suffix is not recognized. Supported formats are ser, mrc, rec, ali, st, emd, dm3, dm4.')

    return out


if __name__ == '__main__':
    dPath = Path(r'c:\users\linol\data')
    out0 = read(r'c:\users\linol\data\10_series.emd')
    out0 = read(r'c:\users\linol\data\1435 1.2 Mx STEM HAADF-DF4-DF2-BF.emd')
    out0 = read(r'c:\users\linol\data\5_Te_15x83_ss=3nm_CL=245_alpha=p48_p06sec_no beamstop_bin4_300kV.dm4')
    fPath = Path('AgNWweld_tomo2_1wire_115kx_160mmCL.mrc')
    out0 = read(dPath / fPath)
    print(out0.keys())