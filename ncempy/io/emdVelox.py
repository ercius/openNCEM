''' Provides an interface to Velox EMD datasets. Not to be confused with
Berkeley EMD data sets (see emd.py) instead.

The reader for EMD Berkeley and Velox files will be combined in the near
future once they are fully tested separately.

Currently limited to only images. 

TODO: Implement h5py_cache to open files with many images
TODO: Test for list_data out of range
'''

import datetime
import json

import numpy as np
import h5py

#import h5py_cache

class fileEMDVelox:
    '''Class to represent Velox EMD files
    
    Parameters:
        filename (str):    Name of the EMD file.
    
    Example
    '''
    
    def __init__(self, filename):
        '''Init opening the file and finding all data groups.
        
        Parameters:
            filename (str): 
        
        Currently only searches the Images group.
        
        '''
        
        ## necessary declaration in case something goes bad
        self.file_hdl = None
        
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string!')

        # try opening the file
        try:
            self.file_hdl = h5py.File(filename, 'r')
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
            out += 'Dataset #{} from detector: {}\n'.format(ii,md['detector'])
        out += 'pixel size = ({0[0]:0.4f}, {0[1]:0.4f}) nm'.format(md['pixelSize'])
        return out
    
    def _find_groups(self):
        '''Find all groups that contain data
        
        Note:
            This currently only finds images
        '''
        try:
            #Get groups and images and metadata
            f1Im = self.file_hdl['Data/Image']
            #Get all of the groups in the Image group
            self.list_data = list(self.file_hdl['Data/Image'].values())
        except:
            self.list_data = []
            raise
    
    def get_dataset(self, group):
        '''Get the data from a group and the associated metadata.
        
        '''
        # check input
        if not isinstance(group, h5py._hl.group.Group):
            raise TypeError('group needs to refer to a valid HDF5 group!')

        data = group['Data'] #the full data set
        metaData = self.parseMetaData(group)
        return (data,metaData)
    
    def parseMetaData(self,group):
        '''Parse metadata in a data group.
        
        Finds the pixelSize and detector name.
        
        '''
        md = {}
        tempMetaData = group['Metadata'][:,0] #get the metadata
        validMetaDataIndex = np.where(tempMetaData > 0) #find valid metadata
        metaData = tempMetaData[validMetaDataIndex].tostring() #change to string
        self.metaDataJSON = json.loads(metaData.decode('utf-8','ignore')) #load UTF-8 string as JSON and output dict
        detectorName = self.metaDataJSON['BinaryResult']['Detector']
        pixelSizeX = float(self.metaDataJSON['BinaryResult']['PixelSize']['width'])*1e9 #change to nm
        pixelSizeY = float(self.metaDataJSON['BinaryResult']['PixelSize']['height'])*1e9 #change to nm
        md['pixelSize'] = (pixelSizeX,pixelSizeY)
        md['detector'] = detectorName
        
        return md
    