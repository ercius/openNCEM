"""
This module provides an interface to Dectris Arina data sets
"""

from pathlib import Path
import h5py
import numpy as np
import hdf5plugin

class fileDECTRIS:
    """ Class to represent Dectris Arina data sets

        Attributes
        ----------
        raw_shape : list
            The shape of the raw data. This is three-dimensional: [num_frames, frameY, frameX].
        data_shape : list
            The four-dimensional shape of the dataset. By default, the 
            scanned region is square.
        file_hdl : h5py.File
            The h5py file handle which provides direct access to the underlying hdf5 file structure.
        data_type : numpy.dtype
            The data type of the values in the data set.
        """
    def __init__(self, filename, bad_pixels=None, verbose=False):
        """ Initialize a data set by opening the master file and determining the file size

            Parameters
            ----------
            filename : str or pathlib.Path or file object
                The HDF5 master file to open.
            verbose : bool, default False
                If True, prints out debugging information
        """
        
        self._verbose = verbose
        self.raw_shape = [0, 0, 0] # shape of data on disk
        self.data_shape = [0, 0, 0, 0] # the shape of the final 4D dataset
        self.file_hdl = None
        self.data_dtype = None
        self.bad_pixel_value = bad_pixels
        
        # Pixels to remove automatically
        # self.bad_pixels = ((49, 75), (93,118), (95,119), (108, 57)) # NCEM bad pixels
        
        if hasattr(filename, 'read'):
            try:
                self.file_path = Path(filename.name)
                self.file_name = self.file_path.name
            except AttributeError:
                self.file_path = None
                self.file_name = None
        else:
            # check filename type, change to pathlib.Path
            if isinstance(filename, str):
                filename = Path(filename)
            elif isinstance(filename, Path):
                pass
            else:
                raise TypeError('Filename is supposed to be a string or pathlib.Path or file object')
            self.file_path = Path(filename)
            self.file_name = self.file_path.name

        # Try opening the file
        try:
            self.file_hdl = h5py.File(filename, 'r')
            assert self.file_hdl['/entry/data']
        except:
            print('Error opening file: "{}"'.format(filename))
            raise

        # if this is a HDF5 file
        if self.file_hdl:
            # Find the initial shape of the data set
            for v in self.file_hdl['/entry/data'].values():
                self.raw_shape[0] = self.raw_shape[0] + v.shape[0]
                self.raw_shape[1] = v.shape[1]
                self.raw_shape[2] = v.shape[2]
                self.data_dtype = v.dtype

    def __del__(self):
        """ Destructor for EMD file object.

        """
        # close the file
        # if(not self.file_hdl.closed):
        self.file_hdl.close()

    def __enter__(self):
        """Implement python's with statement for context managers.

        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Implement python's with statement fr context managers.
        and close the file via __del__()
        """
        self.__del__()
        return None

    def getDataset(self, remove_bad_pixels=False, assume_shape=None):
        """ Read the data from the HDF5 files

        Parameters
        ----------
        remove_bad_pixels : bool, default False
            If True, _remove_bad_pixels function is called after the data is loaded.
        assume_shape : tuple, optional
            If this is set, then this tuple is used as the scanning shape overriding 
            the assumption of a square real space scanning grid
        """
        # Pre allocate space
        data = np.zeros(self.raw_shape, dtype=self.data_dtype)
        # Read in the data in all linked files
        ii = 0
        for v in self.file_hdl['/entry/data'].values():
            data[ii:ii+v.shape[0]] = v[:]
            ii += v.shape[0]

        if assume_shape:
            self.data_shape = (assume_shape[0], assume_shape[1],
                               data.shape[1], data.shape[2])
        else:
            # Reshape assuming square
            shape_square = int((data.shape[0])**0.5)
            assert data.shape[0] == shape_square**2
            self.data_shape = (shape_square, shape_square,
                               data.shape[1], data.shape[2])
        data = data.reshape(self.data_shape)
        if remove_bad_pixels:
            self._remove_bad_pixels()

        data_out = {}
        data_out['data'] = data
        return data_out

    def getMetadata(self):
        """ The dectris Arina files sometimes output an extra file with 
        metadata in it. This checks for that file and reads the meta data
        if if exists. The units are assumed to be nanometers.

        Returns
        -------
        : dict
            Meta data as a dictionary
        
        """

        filename_parts = self.file_path.stem.split('_')
        metadata_file_path = self.file_path.parent / Path('_'.join(filename_parts[0:-1])).with_suffix('.h5')
        if metadata_file_path.exists():
            try:
                metadata = {}
                with h5py.File(metadata_file_path, 'r') as f0:
                    for k,v in f0["STEM Metadata"].attrs.items():
                        metadata[k] = v
                
                pixel_size0 = metadata["Pixel Size"] # convert to ncempy standard
                metadata['pixelSize'] = (pixel_size0, pixel_size0)
                metadata['pixelUnit'] = ('n_m', 'n_m')
                
                return metadata
            except:
                raise
        
    def remove_bad_pixels(self, data, value=0, bad_pixels=None):
        """ Some pixels are known to be very high or very low. This function will replace the 
        pixel values.

        Parameters
        ----------
        data : numpy.ndarray
            The 4D-STEM data set
        value : int or float
            The value to replace the bad pixels by.
        bad_pixels : numpy.ndarray
            A m by 2 ndarray where m is the number of bad pixels and the locations
            are specified in order for frame axis 2 and 3.
        
        """        
        if bad_pixels:
            self.bad_pixels = bad_pixels
        for bad in self.bad_pixels:
            data[:, :, bad[0], bad[1]] = value

def dectrisReader(file_name):
    if isinstance(file_name, str):
        file_name = Path(file_name)

    with fileDECTRIS(file_name) as f1:  # open the file and init the class
        im1 = f1.getDataset()  # read in the dataset
        md = f1.getMetadata()
        if md:
            extra_metadata = {'pixelSize': md['pixelSize'], 'pixelUnit':md['pixelUnit'], 'filename': f1.file_name}
            im1.update(extra_metadata)
    return im1
        