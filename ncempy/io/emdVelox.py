''' Provides an interface to Velox EMD datasets. Not to be confused with
Berkeley EMD data sets (see emd.py) instead.

The reader for EMD Berkeley and Velox files will be combined in the near
future once they are fully tested separately.

Currently limited to only images. This file can not load spectra.
'''

import json

import numpy as np
import h5py

import h5py_cache

class fileEMDVelox:
    '''Class to represent Velox EMD files. It uses the h5py_cache module
    to increase the default cache size from 1MB to 10MB. This significantly
    improves file reading for EMDVelox files which are written with Fortran-
    style ordering and a poor choice of chunking.
    
    Parameters
    ----------
        filename: str
            Name of the EMD file.
    
    Attributes
    ----------
        list_data: list
            A list containing each h5py data group.
        
        file_hdl: h5py.File
            The File handle from h5py.File.
        
        metaDataJSON: dict
            The metadata for the most recently loaded data set.
        
    
    Example
    -------
        Open an EMD Velox file of 1 image.
        
        >>> import ncempy.io as nio
        >>> emd1 = nio.emdVelox.fileEMDVelox('1435 1.2 Mx STEM HAADF-DF4-DF2-BF.emd')
        >>> print(emd1) # print information about the file
        >>> im0, metadata0 = emd1.get_dataset(emd1.list_data[0])
        >>> plt.imshow(im0, extent = (0,md1['pixelSize'][0],0,md1['pixelSize'][1]))
    '''
    
    def __init__(self, filename):
        '''Init opening the file and finding all data groups. Currently only
        searches the /Data/Images group.
        
        Parameters
        ----------
            filename:
                The file path to load as a string.
        
        '''
        
        ## necessary declaration in case something goes wrong
        self.file_hdl = None
        self.metaDataJSON = None
        self.list_data = None
        
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string!')

        # try opening the file
        try:
            #self.file_hdl = h5py.File(filename, 'r')
            self.file_hdl = h5py_cache.File(filename, 'r', chunk_cache_mem_size=10*1024**2)
        except:
            print('Error opening file for readonly: "{}"'.format(filename))
            raise
            
        self._find_groups()
        
    def __del__(self):
        '''Destructor for EMD file object.
        
        Closes the h5py file.
        
        '''
        # close the file
        self.file_hdl.close()

    def __enter__(self):
        '''Implement python's with staement
        
        '''
        return self
        
    def __exit__(self,type,value,traceback):
        '''Implement python's with statment
        and close the file via __del__()
        
        '''
        self.__del__()
        return None
    
    def __str__(self):
        ''' Print out the detectors used to take the data and
        the pixel size to help with telling users about the data in the file.
        
        '''
        out = 'EMD file contains {} data sets\n'.format(len(self.list_data))
        for ii, group in enumerate(self.list_data):
            md = self.parseMetaData(group)
            out += 'Dataset #{} from detector: {}\n'.format(ii,md['detectorName'])
        out += 'pixel size = ({0[0]:0.4f}, {0[1]:0.4f}) nm'.format(md['pixelSize'])
        return out
    
    def _find_groups(self):
        '''Find all groups that contain data.
        
        Note
        ----
            This currently only finds images.
        '''
        try:
            # Get all of the groups in the Image group
            self.list_data = list(self.file_hdl['Data/Image'].values())
        except:
            self.list_data = []
            raise
    
    def get_dataset(self, group):
        '''Get the data from a group and the associated metadata.
        
        Parameters
        ----------
            group: HDF5 dataset or int
                The link to the HDF5 dataset in the file or an integer for the 
                number of the dataset. The list of datasets is held in the
                list_data attribute populated on class init.
        '''
        # check input
        try:
            if type(group) is int:
                group = self.list_data[group]
        except IndexError:
            raise IndexError('EMDVelox group #{} does not exist.'.format(group))
        
        if not isinstance(group, h5py._hl.group.Group):
            raise TypeError('group needs to refer to a valid HDF5 group!')

        data = np.squeeze(group['Data'][:]) #load the full data set
        metaData = self.parseMetaData(group)
        return (data,metaData)
    
    def parseMetaData(self,group):
        '''Parse metadata in a data group. Determines the pixelSize and 
        detector name. The EMDVelox data sets have extensive metadata
        stored as a JSON type string.
        
        Parameters
        ----------
            group: h5py group
                The h5py group to load the metadata from. The string is loaded
                and parsed by the json module into a dictionary
        
        Returns
        -------
            md: dict
                The JSON information returned as a python dictionary.
        
        '''
        md = {}
        tempMetaData = group['Metadata'][:,0]
        # Reduce to valid metadata
        validMetaDataIndex = np.where(tempMetaData > 0) 
        metaData = tempMetaData[validMetaDataIndex].tostring()
        # Interpret as UTF-8 encoded characters and load as JSON
        self.metaDataJSON = json.loads(metaData.decode('utf-8','ignore'))
        # Pull out basic meta data about the images
        detectorName = self.metaDataJSON['BinaryResult']['Detector']
        pixelSizeX = float(self.metaDataJSON['BinaryResult']['PixelSize']['width'])*1e9 #change to nm
        pixelSizeY = float(self.metaDataJSON['BinaryResult']['PixelSize']['height'])*1e9 #change to nm
        # Construct meta data dictionary
        md['pixelSize'] = (pixelSizeX,pixelSizeY)
        md['pixelSizeUnit'] = ('nm','nm')
        md['detectorName'] = detectorName
        
        return md
    