"""
This module provides an interface to SMV files which often end in the IMG format.
These are commonly used in microED.

See https://wiki.uni-konstanz.de/ccp4/index.php/SMV_file_format for more details.

Note
----
General users:
    Use the simplified smv.smvReader() function to load the data and meta
    data as a python dictionary.

Advanced users and developers:
    Access the file internals through the smv.fileSMV() class.

"""

from pathlib import Path
import datetime

import numpy as np

class fileSMV:
    """Class to represent SMV files.

    Attributes
    ----------
    file_name : str
        The name of the file
    file_path : pathlib.Path
        A pathlib.Path object for the open file
    camera_length=110, lamda=0.0197, pixel_size=0.01, beam_center=None, binned_by=1
    """
    
    def __init__(self, filename, verbose=False):
        """ Initialize opening the file
        
        Parameters
        ----------
        filename : str or pathlib.Path or file object
            String pointing to the filesystem location of the file.

        verbose : bool, optional, default False
            If True, debug information is printed.
            
        """
        
        self._verbose = verbose
        self._expected_keys = ('HEADER_BYTES', 'DIM', 'BYTE_ORDER', 'TYPE', 'SIZE1', 'SIZE2', 'PIXEL_SIZE', 
                         'WAVELENGTH', 'DISTANCE', 'PHI', 'BEAM_CENTER_X', 'BEAM_CENTER_Y', 'BIN', 
                         'DATE', 'DETECTOR_SN', 'OSC_RANGE', 'OSC_START', 'IMAGE_PEDESTAL', 'TIME', 
                         'TWOTHETA')
        self.header_info = {}
        self.num_header_bytes = 512
        
        if hasattr(filename, 'read'):
            self.fid = filename
            try:
                    self.file_name = self.fid.name
            except AttributeError:
                self.file_name = None
        else:
            # check filename type. Prefer pathlib.Path
            if isinstance(filename, str):
                filename = Path(filename)
            elif isinstance(filename, Path):
                pass
            else:
                raise TypeError('Filename is supposed to be a string or pathlib.Path')
        
        self.file_path = filename
        self.file_name = self.file_path.name
        
        try:
            self.fid = open(self.file_path, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(self.file_path))
            raise
        except:
            raise
        
        if not self._validate():
            raise IOError('Not an SMV file: {}'.format(self.file_path))
        
    def __del__(self):
        """Destructor which also closes the file

        """
        if not self.fid.closed:
            if self._v:
                print('Closing input file: {}'.format(self.file_path))
            self.fid.close()

    def __enter__(self):
        """Implement python's with statement

        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Implement python's with statement
        and close the file via __del__()
        """
        self.__del__()
        return None
    
    def _validate(self):
        first_line = self.read(15).decode('UTF-8')
        if first_line = '{\nHEADER_BYTES=':
            bytes_str = self.fid.readline()
            self.num_header_bytes = int(bytes_str.strip().strip(';'))
            return True
        else:
            return False
    
    def readHeader(self):
        """Read the header information and conver to numbers or strings."""
        
        self.fid.seek(0, 0)
        head = f0.read(self.num_header_bytes).decode('UTF-8').split('\n')
        for line in head:
            if '=' in line:
                key, val = line.split('=')
                val = val.strip(';')
                if key in self._expected_keys:
                    try:
                        self.header_info[key] = float(val)
                        try:
                            self.header_info[key] = int(val)
                        except:
                            pass # it is a float
                    except:
                        self.header_info[key] = val # not a number
    
    def parseHeader():
        """Parse the header dictionary for relelvant information to read the data in the file."""
        for key, val in self.header_info:
            if key == ''
    
    def getDataset(self):
        self.readHeader()
        self.parseHeader()
        self.getDataset()
        
def smvWriter(out_path, dp, camera_length=110, lamda=0.0197, pixel_size=0.01, beam_center=None, binned_by=1):
    """ Write out diffraction as SMV formatted file
    Header is 512 bytes of zeros and then filled with ASCII
    
    only little endian is supported.
    
    only uint16 is supported.
    
    camera length in mm
    lamda = 0.0197  # angstroms
    pixel_size = 0.01  # physical detector pixel size in mm
    beam_center in column, row format in mm
    """
    print('New write_smv version!!')
    
    if dp.dtype != np.uint16:
        raise TypeError("Only uint16 data type is supported.")
    dtype = 'unsigned_short'
    
    if not beam_center:
        beam_center = [ii / 2 * pixel_size for ii in dp.shape]
    
    # make sure binned_by is an integer
    binned_by = int(binned_by)
    
    # Write 512 bytes of zeros
    with open(out_path, 'wb') as f0:
        f0.write(np.zeros(512, dtype=np.uint8))
    # Write the header over the zeros as needed
    with open(out_path, 'r+') as f0:
        f0.write("{\nHEADER_BYTES=512;\n")
        f0.write("DIM=2;\n")
        f0.write("BYTE_ORDER=little_endian;\n")
        f0.write(f"TYPE={dtype};\n")
        f0.write(f"SIZE1={dp.shape[1]};\n")  # size 1 is columns
        f0.write(f"SIZE2={dp.shape[0]};\n")  # size 2 is rows
        f0.write(f"PIXEL_SIZE={pixel_size};\n")  # physical pixel size in micron
        f0.write(f"WAVELENGTH={lamda};\n")  # wavelength
        f0.write(f"DISTANCE={int(camera_length)};\n")
        f0.write("PHI=0.0;\n")
        f0.write(f"BEAM_CENTER_X={beam_center[1]};\n")
        f0.write(f"BEAM_CENTER_Y={beam_center[0]};\n")
        f0.write(f"BIN={binned_by}x{binned_by};\n")
        f0.write("DATE=Fri Dec 31 23:59:59 1999;\n")
        f0.write("DETECTOR_SN=unknown;\n")
        f0.write("OSC_RANGE=1.0;\n")
        f0.write("OSC_START=0;\n")
        f0.write("IMAGE_PEDESTAL=0;\n")
        f0.write("TIME=10.0;\n")
        f0.write("TWOTHETA=0;\n")
        f0.write("}\n")
    # Append the binary image data at the end of the header
    with open(out_path, 'rb+') as f0:
        f0.seek(512, 0)
        f0.write(dp)