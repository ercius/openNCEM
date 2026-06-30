"""
This module provides an interface to the FEI / Thermo Fischer ESVision EMI files. 
It is only meant to open files that do not have a corresponding .ser file.

Basic EMI reader which can read the first image in most emi files.
Currently readable image data types are uint16, uint32, float

"""

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
        self.file_path = filename
        self.file_name = self.file_path.name

        try:
            self.fid = open(self.file_path, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(self.file_path))
            raise
        except:
            raise

    def _read_text(self, fid):
        aa = np.fromfile(fid, dtype=self._text_dtype, count=1)
        if aa['size'] > 0:
            bin = np.fromfile(fid, dtype='<u1', count=aa['size'][0])
            text = ''.join([chr(item) for item in bin])
        else:
            text = ''
        return text

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
    
    def parse_file(self):
        # Read in the full file
        full_file = np.fromfile(self.file_path, dtype='<u1')
        # Find bytes that indicate a certain ASCII character: `
        obj_loc = np.where(full_file == 96)[0]
        
        with open(self.file_path, 'rb') as f0:
            for loc in obj_loc[0:]:
                f0.seek(loc, 0)
                cur_text = self._read_text(f0)
                self.image_name.append(cur_text)
                if self._verbose:
                    print(cur_text)
                second_field = np.fromfile(f0, dtype='<u2', count=1)
                if (second_field == 112) or (second_field == 17184):
                    obj_info = np.fromfile(f0, count=2, dtype='<u2')
                    if obj_info[0] == 1042:
                        # this is an image. Read header and continue
                        image_info = np.fromfile(f0, count=10, dtype='<u2')
                        f0.seek(-8, 1);
                        self.image_size.append(np.fromfile(f0,count=2,dtype='<u4'))
                        if image_info[3] == 8710:
                            self.data_type.append('<u2')
                        elif image_info[3] == 8714:
                            self.data_type.append('<u4')
                        elif image_info[3] == 514:
                            self.data_type.append('<u4')
                        elif image_info[3] == 8716:
                            self.data_type.append('f32')
                        else:
                            print('Unknown data type: {}'.format(image_info[3]))
                            print('for object named: {}'.format(cur_text))
                            return
                        self.image_locations.append(self.fid.tell())
                        
                f0.seek(-2, 1) # roll back the pointer 2 bytes
    
    def getDataset(self, index=0):
        """Read the data from the file
        
        Paremeters
        ----------
        index : int, optional
        The index of the image to load.

        Returns
        -------
        : dict
        A dictionary contraining the data with the key 'data'
        
        """
        self.fid.seek(self.image_locations[index])
        image_size = self.image_size[index]
        dtype = self.data_type[index]
        image = np.fromfile(self.fid, count=image_size[0]*image_size[1], dtype=self.data_type)
        print('Read image named: {}'.format(self.image_name[index]))

        return {'data': image.reshape(image_size)}

def emiReader(fname):
    """Old original function to read EMI files. Use the fileEMI class instead."""
    full_file = np.fromfile(fname,dtype='<u1')
    
    obj_loc = np.where(full_file == 96)[0] # character `
    
    with open(fname,'rb') as f0:
        for loc in obj_loc[0:]:
            f0.seek(loc, 0)
            cur_text = read_text(f0)
            if v:
                print(cur_text)
            second_field = np.fromfile(f0,dtype='<u2',count=1)
            if (second_field == 112) or (second_field == 17184):
                obj_info = np.fromfile(f0,count=2,dtype='<u2')
                if obj_info[0] == 1042:
                    #this is the data. Read header and quit
                    image_info = np.fromfile(f0,count=10,dtype='<u2')
                    f0.seek(-8,1);
                    im_size = np.fromfile(f0,count=2,dtype='<u4')
                    if image_info[3] == 8710:
                        numType = '<u2'
                    elif image_info[3] == 8714:
                        numType = '<u4'
                    elif image_info[3] == 514:
                        # Flucam data set
                        numType = '<u4'
                    elif image_info[3] == 8716:
                        numType = 'f32'
                    else:
                        print('Unknown data type: {}'.format(image_info[3]))
                        print('for object named: {}'.format(cur_text))
                        return
                    
                    image = np.fromfile(f0,count=im_size[0]*im_size[1],dtype=numType)
                    print('Only read first image named: {}'.format(cur_text))
                    print(image)
                    return image.reshape(im_size)
            f0.seek(-2, 1) # roll back the pointer 2 bytes
    
if __name__ == '__main__':
    im = emiReader('c:/users/linol/data/with sample CL 170.emi')
