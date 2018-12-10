''' Provides an interface to Velox EMD datasets. Not to be confused with
Berkeley EMD data sets (see emd.py) instead.

The reader for EMD Berkeley and Velox files will be combined in the near
future once they are fully tested separately.

Currently limited to only images. 
'''

import numpy as np
import h5py
import datetime


class fileEMDVelox:
    '''Class to represent Velox EMD files
    
    Parameters:
        filename (str):    Name of the EMD file.
        readonly (bool):    Set to open in read only mode.
    
    '''
    
    def __init__(self, filename, readonly=False):
        '''Init opening/creating the file.
        
        '''
        
        ## necessary declarations in case something goes bad
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
            
        self.find_groups()
        self.parseMetaData()
        
        '''
        #Plot the images
        if len(dsetGroups) > 1:
            fg1, ax1 = plt.subplots(nrows=round(len(dsetGroups)/2+0.1),ncols=2,figsize=(8,12)) #two columns
            for ii,image in enumerate(dsetGroups):
                fovX = pixelSizeX[ii]*image['Data'].shape[0]
                fovY = pixelSizeY[ii]*image['Data'].shape[1]
                ax1.ravel()[ii].imshow(image['Data'][:,:,0],extent=[0,fovX,0,fovY],origin='lower')
                ax1.ravel()[ii].set(xlabel='X (nm)',ylabel='Y (nm)',title=detectorName[ii])
        else:
            fg1, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(8,12)) #two columns
            fovX = pixelSizeX[ii]*image['Data'].shape[0]
            fovY = pixelSizeY[ii]*image['Data'].shape[1]
            ax1.imshow(image['Data'][:,:,0],extent=[0,fovX,0,fovY],origin='lower')
            ax1.set(xlabel='X (nm)',ylabel='Y (nm)',title=detectorName[ii])
        fg1.tight_layout()
        '''
        
    def __del__(self):
        '''Destructor for EMD file object. 
        
        '''
        # close the file
        #if(not self.file_hdl.closed):
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
    
    def find_groups(self):
        '''Find all groups that contain data
        
        Note:
            This currently only finds images
        '''
        try:
            #Get groups and images and metadata
            f1Im = file_hdl['Data/Image']
            #Get all of the groups in the Image group
            dsetGroups = list(file_hdl['Data/Image'].values())
        except:
            dsetGroups = []
            raise
        return dsetGroups
    
    def parseMetaData(self):
        '''Parse metadata in the file.
        
        FInds the pixelSize and detector name. All metadata is a string saved as a
        class variable called metaData.
        
        '''
        
        self.pixelSizeX = []
        self.pixelSizeY = []
        self.detectorName = []
        for image in dsetGroups:
            tempMetaData = image['Metadata'][:,0] #get the metadata
            validMetaDataIndex = np.where(tempMetaData > 0) #find valid metadata
            self.metaData = tempMetaData[validMetaDataIndex].tostring() #change to string
            self.metaDataJSON = json.loads(self.metaData.decode('utf-8','ignore')) #load UTF-8 string as JSON and output dict
            self.detectorName.append(self.metaDataJSON['BinaryResult']['Detector'])
            self.pixelSizeX.append(float(self.metaDataJSON['BinaryResult']['PixelSize']['width'])*1e9) #change to nm
            self.pixelSizeY.append(float(self.metaDataJSON['BinaryResult']['PixelSize']['height'])*1e9) #change to nm
    
    def get_dataset(self, groupNum):
        '''Get the data from a group and the associated metadata.
        
        '''
        
        data = self.dsetGroups[groupNum]['Data'] #the full data set
        metaData = {}
        metadata['pixelSize' = (self.pixelSizeX[groupNum]
        metaData['detector'] = self.detectorName[groupNum]
        return (data,metaData)
        