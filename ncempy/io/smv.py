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
    fid : file
        The file handle to the opened MRC file.
    dataType : np.dtype
        The numpy dtype of the data.
    dataSize : np.ndarray
        The number of pixels along each dimension. Corresponds to the shape attribute of a np.ndarray
    num_header_bytes : 
        The number of bytes in the header. Usually 512.
    header_info : dict
        A dictionary containing the header meta data.

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
        self._data_types = {'unsigned_short': np.uint16} # convert SMV types with numpy dtypes
        self.header_info = {}
        self.num_header_bytes = None
        self.dataType = None
        self.dataSize = [0, 0]
        self._v = verbose
        
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
        
        self.readHeader()
        self.parseHeader()
        
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
        # first_line = self.fid.read(15).decode('UTF-8')
        first_line = self.fid.readline() # should be {
        second_line = self.fid.readline().decode('UTF-8') # should contain HEADER_BYTES
        if 'HEADER_BYTES=' in second_line:
            self.num_header_bytes = int(second_line.split('=')[1].strip().strip(';'))
            return True
        else:
            return False
    
    def readHeader(self):
        """Read the header information and conver to numbers or strings."""
        
        self.fid.seek(0, 0)
        head = self.fid.read(self.num_header_bytes).decode('UTF-8').splitlines()
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
    
    def parseHeader(self):
        """Parse the header dictionary for relelvant information to read the data in the file."""
        for key, val in self.header_info.items():
            if key == 'SIZE1':
                self.dataSize[1] = val # column
            elif key == 'SIZE2':
                self.dataSize[0] = val # row
            elif key == 'TYPE':
                try:
                    self.dataType = self._data_types[val]
                except KeyError:
                    raise(f'File data type not supported: {val}')
    
    def getDataset(self):
        self.readHeader()
        self.parseHeader()
        
        self.fid.seek(self.num_header_bytes, 0)
        data = np.fromfile(self.fid, count=self.dataSize[0] * self.dataSize[1], dtype=self.dataType)
        data = data.reshape(self.dataSize)
        data_out = {}
        data_out['data'] = data
        return data_out
    
def smvWriter(out_path, dp, camera_length=110, lamda=0.0197, pixel_size=0.01, beam_center=None, binned_by=1, newline=None):
    """ Write out data as a SMV (.img) formatted file
    Header is 512 bytes of zeros and then filled with ASCII.
    
    Note:
    - only little endian is supported
    - only uint16 is supported
    - ony 2D data is supported
    - some other meta data (PHI, DATE, etc.) is populated with hard coded values
    
    Parameters
    ----------
    camera_length : float
        The calibrated camera length (not the label) in mm. Default is 110 mm.
    lamda : float
        The wavelength of the radiation in Ansgroms. Default is 0.0197 for 300 kV electrons
    pixel_size : float
        Physical detector pixel size in mm. Default is 0.01 mm (10 microns)
    beam_center : tuple
        The location of the center beam in column, row format in mm (not pixels!)
    binned_by : int
        The binning applied to the original data. This is necessary for proper
        calibrations of detector distances and beam center. Default is 1.
    newline : str (optional)
        Allow the user to specify the newline character. For data written on Windows computers
        some processing programs in Linux are not able to load SMV files with Windows 
        carriage return and newline characters. Use '\n' on Windows machines to enforce Linux
        line endings. The default `None` will use the system default.
    """
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
    # The newline character is system default unless otherwise specified
    with open(out_path, 'r+', newline=newline) as f0:
        f0.write("{\nHEADER_BYTES=512;\n")
        f0.write("DIM=2;\n")
        f0.write("BYTE_ORDER=little_endian;\n")
        f0.write(f"TYPE={dtype};\n")
        f0.write(f"SIZE1={dp.shape[1]};\n")  # size 1 is columns
        f0.write(f"SIZE2={dp.shape[0]};\n")  # size 2 is rows
        f0.write(f"PIXEL_SIZE={pixel_size};\n")  # physical pixel size in mm
        f0.write(f"WAVELENGTH={lamda};\n")  # wavelength in Angstroms
        f0.write(f"DISTANCE={int(camera_length)};\n") # in mm
        f0.write("PHI=0.0;\n")
        f0.write(f"BEAM_CENTER_X={beam_center[1]};\n") # in mm (not pixels!)
        f0.write(f"BEAM_CENTER_Y={beam_center[0]};\n") 
        f0.write(f"BIN={binned_by}x{binned_by};\n")
        f0.write("DATE=Fri Dec 31 23:59:59 1999;\n")
        f0.write("DETECTOR_SN=unknown;\n")
        f0.write("OSC_RANGE=1.0;\n")
        f0.write("OSC_START=0;\n")
        f0.write("IMAGE_PEDESTAL=0;\n")
        f0.write("TIME=1.0;\n")
        f0.write("TWOTHETA=0;\n")
        f0.write("}\n")
    # Append the binary image data at the end of the header
    with open(out_path, 'rb+') as f0:
        f0.seek(512, 0)
        f0.write(dp)

        
def smvReader(file_name, verbose=False):
    """ A simple function to read open a SMV, parse the header, and read the
    data and meta data.

    Parameters
    ----------
        file_name : str or pathlib.Path
            The path to the file to load.

    Returns
    -------
        out : dict
            A dictionary containing the data and interesting metadata.
    
    Note
    ----
    The returned dictionary has an entry called `pixelSize`. This is the calibrated distance
    in angstroms. It is not the physical size of a detector pixel. 

    Example
    -------
        Simply read in all data from disk into memory. This assumes the dataset is 3 dimensional:
        >> from ncempy.io.smv import smvReader
        >> import matplotlib.pyplot as plt
        >> mrc1 = smvReader.('filename.mrc')
        >> plt.imshow(mrc1['data'][0, :, :]) #show the first image in the data set
    """
    if isinstance(file_name, str):
        file_name = Path(file_name)

    with fileSMV(file_name) as f1:  # open the file and init the class
        im1 = f1.getDataset()  # read in the dataset

        # Calculate the pixel size in inverse angstroms according to the geometry in the header
        BIN = [int(ii) for ii in f1.header_info['BIN'].split('x')]
        alpha = (BIN[0] * f1.header_info['PIXEL_SIZE']) / f1.header_info['DISTANCE'] # angle across 1 pixel
        dp_pixel_distance = alpha / f1.header_info['WAVELENGTH'] * 1e-10 # divide by wavelength to get distance in Angstroms
        pixelSize = (dp_pixel_distance, dp_pixel_distance)
        extra_metadata = {'pixelSize': pixelSize, 'pixelUnit':'A', 'filename': f1.file_name}
    im1.update(extra_metadata)
    return im1
