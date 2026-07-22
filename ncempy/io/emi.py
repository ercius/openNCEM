"""
This module provides an interface to the FEI / Thermo Fischer ESVision EMI files. 
It is only meant to open files that do not have a corresponding .ser file. It can 
read the first image in some emi files.

Currently readable image data types are uint16, uint32, float32

The EMI file format is not documented at all. This reader is based on 
reverse engineering of a few sample files from only a couple of microscopes. 
It is likely that there are other variations of the EMI file format that this reader will not be able to read.

"""

import mmap
from pathlib import Path
import numpy as np

class fileEMI:
    """Class to represent simple EMI files.

    Attributes
    ----------
    file_name : str
        The name of the file
    file_path : pathlib.Path
        A pathlib.Path object for the open file
    fid : file
        The file handle to the opened file.
    data_type : list of np.dtype
        The numpy dtype of the data.
    image_size : list of 2-tuples
        The size of the images in pixels.
    image_locations : list of int
        The file bytes locations of the images in the file.
    image_name : list of str
        The names of the images in the file.
    """
    
    _text_dtype = np.dtype([('mark','<u2'),('unknown','<u2'),('size','<u4')])

    def __init__(self, file_name, verbose=False):
        
        self.data_type = []
        self.image_size = []
        self.image_locations = []
        self.image_name = []
        self._verbose = verbose

        if hasattr(file_name, 'read'):
            self.fid = file_name
            try:
                self.file_name = self.fid.name
            except AttributeError:
                self.file_name = None
        else:
            # check filename type. Prefer pathlib.Path
            if isinstance(file_name, str):
                filename = Path(file_name)
            elif isinstance(file_name, Path):
                pass
            else:
                raise TypeError('Filename is supposed to be a string or pathlib.Path')
        self.file_path = file_name
        self.file_name = self.file_path.name

        try:
            self.fid = open(self.file_path, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(self.file_path))
            raise
        except:
            raise

    def __del__(self):
        """Destructor which also closes the file

        """
        if not self.fid.closed:
            if self._verbose:
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
    
    def parse_file(self):
        """Parse the file to find the image location, data type, and size. This
        is reverse engineered and only works for very simple files.
        """
        # Map the entire file (0 = entire file)
        with mmap.mmap(self.fid.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            while True:
                pos = mm.find(bytes([96])) # this indicates the start of a data object
                if pos != -1:
                    mm.seek(pos)  # Set mmap pointer to byte location
                    raw_bytes = mm.read(8)
                    text_mark = int.from_bytes(raw_bytes[0:2], byteorder='little', signed=False)
                    _ = int.from_bytes(raw_bytes[2:4], byteorder='little', signed=False)
                    text_size = int.from_bytes(raw_bytes[4:8], byteorder='little', signed=False)
                    if text_size > 0:
                        cur_text = mm.read(text_size).decode('ascii')
                    self.image_name.append(cur_text)
                    if self._verbose:
                        print(cur_text)
                    raw_bytes = mm.read(2)
                    second_field = int.from_bytes(raw_bytes, byteorder='little', signed=False)
                    if (second_field == 112) or (second_field == 17184):
                        raw_bytes = mm.read(4) # read 4 bytes but only convert the first 2
                        obj_info = int.from_bytes(raw_bytes[0:2], byteorder='little', signed=False)
                        if self._verbose:
                            print(f'obj_info byte = {obj_info}')
                        if obj_info == 1042:
                            # this is an image. Read header and continue
                            raw_bytes = mm.read(20)
                            image_dtype = int.from_bytes(raw_bytes[4:6], byteorder='little', signed=False)
                            image_size = (int.from_bytes(raw_bytes[12:17], byteorder='little', signed=False),
                                            int.from_bytes(raw_bytes[16:20], byteorder='little', signed=False))
                            if self._verbose:
                                print(f'image size = {image_size}')
                            self.image_size.append(image_size)
                            if image_dtype == 8710:
                                self.data_type.append('<u2')
                            elif image_dtype == 8714:
                                self.data_type.append('<u4')
                            elif image_dtype == 514:
                                self.data_type.append('<u4')
                            elif image_dtype == 8716:
                                self.data_type.append('<f4')
                            else:
                                print('Unknown data type: {}'.format(image_dtype))
                                print('for object named: {}'.format(cur_text))
                                return
                            if self._verbose:
                                print(f'image location = {mm.tell()}')
                            
                            self.image_locations.append(mm.tell())
                            # We stop at the first image
                            break

    
    def getDataset(self, index=0):
        """Read the data from the file
        
        Paremeters
        ----------
        index : int, optional
        The index of the image to load.

        Returns
        -------
        : dict
        A dictionary containing the data with the key 'data'
        
        """
        self.fid.seek(self.image_locations[index], 0)
        image_size = self.image_size[index]
        dtype = self.data_type[index]
        image = np.fromfile(self.fid, count=image_size[0]*image_size[1], dtype=dtype)
        if self._verbose:
            print('Read image named: {}'.format(self.image_name[index]))

        return {'data':image.reshape(image_size), 'name':self.image_name[index]}

if __name__ == '__main__':
    with fileEMI('../data/with sample CL 170.emi', verbose=False) as f0:
        f0.parse_file()
        data = f0.getDataset(0)
        print(f'data shape: {data["data"].shape}')
