'''
A module to read MRC files in python and numpy.
Written according to MRC specification at http://bio3d.colorado.edu/imod/betaDoc/mrc_format.txt
ALso works with FEI MRC files which include a special header block with experimental information.
written by: Peter Ercius, percius@lbl.gov
'''

import numpy as np

class fileMRC:
    
    def __init__(self, filename, verbose = False):
        '''Init opening the file and reading in the header.
        Read in the data in MRC format and other useful information.
                
        Parameters:
            filename (str): string pointing to the filesystem location of the file.
            verbose (bool): if True, debug information is printed.
        Returns:
            (dict): A dictionary with keys stack, voxelSize, filename, axisOrientations, {FEIinfo}
        
        Note:
            Most users will prefer to use the mrc.mrcReader() function to simply read
            the entire data set into memory with a single command.
        
        '''
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string')
            
        self.filename = filename

        # necessary declarations, if something fails
        self.fid = None
        
        self.dataOut = {} #will hold the data and metadata to output to the user after getDataset() call
        
        #Add a top level variable to indicate verbosee output for debugging
        self.v = verbose
        
        #Open the file and quit if the file does not exist
        try:
            self.fid = open(self.filename,'rb')
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
        
        #Store the original filename
        self.dataOut['filename'] = self.filename
        
        return None
    
    def __del__(self):
        '''Close the file.
        
        '''
        if(not self.fid):
            if self.v:
                print('Closing input file: {}'.format(self.filename))
            self.fid.close()
        return None
    
    def parseHeader(self):
        '''Read the header information which includes data type, data size, data
        shape, and metadata.
        
        Note:
            This header uses Fortran-style ordering. Numpy uses C-style ordering.
            The header is read in and then reversed [::-1] at the end for output
            to the user.
        
        '''
        #Read in the initial header values
        head1 = np.fromfile(self.fid,dtype=np.int32,count=10)
        if self.v:
            print('header1 = {}'.format(head1))
        #Set the number of pixels for each dimension
        self.dataSize = head1[0:3]
        if self.v:
            print('dataSize (fortran ordering) = {}'.format(self.dataSize))
        
        #Set the data type and convert to numpy type
        self.mrcType = head1[3]
        self.dataType = self._getMRCType(self.mrcType)
        if self.v:
            print('dataType = {}'.format(self.dataType))
        
        #Get the grid size
        self.gridSize = head1[7:10]
        if self.v:
            print('mrc defined gridSize = {}'.format(self.gridSize))
            
        #Get the physical volume size (always in Angstroms) (starting at byte #11 in the file).
        head2 = np.fromfile(self.fid,dtype=np.float32,count=6)
        
        self.volumeSize = head2[0:3]
        if self.v:
            print('mrc defined volumeSize = {}'.format(self.volumeSize))
        
        #calculate the voxel size based on volume and grid sizes
        if self.volumeSize.any() ==0 or self.gridSize.any() ==0:
            if self.v:
                print('Detected 0 volume or grid size. Setting voxel size to 1 (Ang).')
            self.voxelSize = np.ones(3) #use 1 as a voxel size if its not set in the file
        else:
            self.voxelSize = (self.volumeSize / np.float32(self.gridSize))
            if self.v:
                print('voxelSize (Ang) = {}'.format(self.voxelSize))

        #Pixel (cell) angles
        self.cellAngles = head2[3:6]
        if self.v:
            print('cellAngles = {}'.format(self.cellAngles))
            
        #Axis orientations. Tells which axes are X,Y,Z
        self.axisOrientations = np.fromfile(self.fid,dtype=np.int32,count=3)
        if self.v:
            print('axisOrientations = {}'.format(self.axisOrientations))
        
        self.Shape = [self.dataSize[x-1] for x in self.axisOrientations]
        if self.v:
            print('data shape = {}'.format(self.Shape))
        
        #Min, max,mean
        self.minMaxMean = np.fromfile(self.fid,dtype=np.int32,count=3)
        
        #Extra information (for FEI MRC file, extra(1) is the size of the FEI information encoded with the file in terms of 4 byte floats)
        self.extra = np.fromfile(self.fid,dtype=np.int32,count=34)
        
        #Numpy uses C-style ordering. The header is written in Fortran-Style ordering. Flip the order of everything useful
        if self.v:
            print('Note: The MRC header is written in Fortran-Style ordering, but Numpy uses C-style ordering. This program will now reverse the order (using [::-1]) of useful metadata: (dataSize, gridSize,volumeSize,voxelSize,cellAngles,axisOrientations)')
        self.dataSize = self.dataSize[::-1]
        self.gridSize = self.gridSize[::-1]
        self.volumeSize = self.volumeSize[::-1]
        self.voxelSize = self.voxelSize[::-1]
        self.cellAngles = self.cellAngles[::-1]
        self.axisOrientations = self.axisOrientations[::-1]
        self.Shape = self.Shape[::-1]
        
        #Move to the end of the normal header
        self.fid.seek(1024)
        
        #Read in the extended header if it exists (for FEI MRC files)
        if self.extra[1] != 0:
            pos1 = self.fid.tell()
            if self.v:
                print('Extra header found. Most likely and FEI-style MRC file.')
                print('Position before reading extra header: ' + str(pos1))
                print('Extra header size = ' +str(self.extra[1]))
            
            '''Read the extra FEI header described as follows:
             1 a_tilt  first Alpha tilt (deg)
             2 b_tilt  first Beta tilt (deg)
             3 x_stage  Stage x position (Unit=m. But if value>1, unit=???m)
             4 y_stage  Stage y position (Unit=m. But if value>1, unit=???m)
             5 z_stage  Stage z position (Unit=m. But if value>1, unit=???m)
             6 x_shift  Image shift x (Unit=m. But if value>1, unit=???m)
             7 y_shift  Image shift y (Unit=m. But if value>1, unit=???m)
             8 defocus  starting Defocus Unit=m. But if value>1, unit=???m)
             9 exp_time Exposure time (s)
             10 mean_int Mean value of image
             11 tilt_axis   Tilt axis (deg)
             12 pixel_size  Pixel size of image (m)
             13 magnification   Magnification used
             14 voltage accelerating voltage
             15 ??
            '''
            FEIinfoValues = np.fromfile(self.fid,dtype=np.float32,count=15)
            self.FEIinfo = {'a_tilt':FEIinfoValues[0],'b_tilt':FEIinfoValues[1],'x_stage':FEIinfoValues[2],'y_stage':FEIinfoValues[3],'z_stage':FEIinfoValues[4],'x_shift':FEIinfoValues[5],'y_shift':FEIinfoValues[6],'defocus':FEIinfoValues[7],'exposure_time':FEIinfoValues[8],'mean':FEIinfoValues[9],'tilt_axis':FEIinfoValues[10],'pixel_size':FEIinfoValues[11],'magnification':FEIinfoValues[12],'voltage':FEIinfoValues[13],'unknown':FEIinfoValues[14]}
            
            self.voxelSize[0] = 1. #set this to 1 but it should be the tilt angles. These can be non-uniform though.
            self.voxelSize[1] = self.FEIinfo['pixel_size']*1e10 #convert [m] to Angstroms as is the standard for MRCs
            self.voxelSize[2] = self.FEIinfo['pixel_size']*1e10
            
            if self.v:
                print('Extended header data')
                for aa,bb in self.FEIinfo.items():
                    print('{} = {}'.format(aa,bb))
            
        self.dataOffset = 1024+self.extra[1] #offset of the data from the start of the file
        
        #Add relevant information (metadata) to the output dictionary
        self.dataOut = {'voxelSize':self.voxelSize,'axisOrientations':self.axisOrientations,'cellAngles':self.cellAngles,'axisOrientations':self.axisOrientations}
        if self.extra[1] != 0:
            self.dataOut['FEIinfo'] = self.FEIinfo
        
        return 1
    
    def getDataset(self):
        '''Read in the full data block and reshape to a matrix
        with C-style ordering.
        
        '''
        self.fid.seek(self.dataOffset,0) #move to the start of the data from the start of the file
        try:
            data1 = np.fromfile(self.fid,dtype=self.dataType,count=np.prod(self.dataSize))
            self.dataOut['data'] = data1.reshape(self.Shape)
        except MemoryError:
            print("Not enough memory to read in the full data set")        
        return self.dataOut
    
    def getSlice(self,num):
        '''Read in a slice of an MRC file. Useful for parsing through a large file without reading
        the entire data set into memory.
        
        Paremeters:
            num (int): Get the requested image.
        
        Returns:
            (ndarray): A 2D slice or a 3D set of slices along the first index
        
        Raises:
            IndexError: If num > the number of slices.
        '''

        # Check num is within the data array size bounds
        if num > (self.dataSize[0]-1):
            raise IndexError('Index {} is out of bounds for array with size {}'.format(num, self.dataSize[0]))
        
        self.fid.seek(self.dataOffset,0) #move to the start of the data by skipping the header
        imSize = self.dataSize[1]*self.dataSize[2] #size of each image in pixels
        byteSize = np.dtype(self.dataType).itemsize*imSize
        self.fid.seek(num*byteSize,1) #skip to the slice requested from the start of the data
        
        data1 = np.fromfile(self.fid,dtype=self.dataType,count=imSize) #read in the requested image
        data1 = data1.reshape((self.Shape[1],self.Shape[2])) #reshape the image
            
        return data1
    
    def _applyAxisOrientations(self,arrayIn):
        ''' This is untested and unused.
        
        '''
        return [arrayIn[x-1] for x in self.axisOrientations]
    
    def _getMRCType(self, dataType):
        '''Return the correct data type according to the official MRC type list:
        0 image : signed 8-bit bytes range -128 to 127
        1 image : 16-bit halfwords
        2 image : 32-bit reals
        3 transform : complex 16-bit integers
        4 transform : complex 32-bit reals
        6 image : unsigned 16-bit range 0 to 65535
        
        Parameters:
            dataType (int): The data type value encoded in an MRC header
            
        Returns:
            {numpy dtype}: The corresponding numpy data type.
        
        '''
        if dataType == 0:
            Type = np.int8
        elif dataType == 1:
            Type = np.int16
        elif dataType == 2:
            Type = np.float32
        elif dataType ==  6:
            Type = np.uint16
        else:
            print("Unsupported data type" + str(dataType)) #complex data types are currently unsupported
        return Type
#end class fileMRC

def mrcReader(fname,verbose=False):
    '''A simple function to read open a MRC, parse the header, and read the full
    data set.
    
    Parameters:
        fname (str): The name of the file to load
        
    Keywords:
        verbose (bool): Enable printing debug messages as the header is parsed.
        
    Returns:
        (dict): A dictionary containing the data and interesting metadata. The data
        is attached to the 'data' key.
    '''
    f1 = fileMRC(fname,verbose) #open the file and init the class
    f1.parseHeader() #parse the header
    im1 = f1.getDataset() #read in the dataset
    del f1 #delete the class and close the file
    return im1 #return the data and metadata as a dictionary

def mrc2raw(fname):
    '''Writes the image data in an MRC file as binary file with the same file
    name and .raw ending. Data type and size are written in the file name.
    No other header information is retained.
    
    Parameters:
        fname (str): The name of the file to convert
    
    '''
    tomo = mrcReader(fname)
    rawName = tomo['filename'].rsplit('.',1)[0] + '_' + str(tomo['stack'].dtype) + '_' + str(tomo['stack'].shape) + '.raw'
    fid = open(rawName,'wb')
    fid.write(tomo['stack']) #write out as C ordered data
    fid.close()
    
def mrc2emd(fname):
    '''Write an MRC file as an HDF5 file in EMD format with same file name and .emd ending.
    Header information is retained as attributes.
    
    TODO: Update this to use ncempy.emd class
    '''
    import h5py
    
    #Read in the MRC data and reshape to C-style ordering
    tomo = mrcReader(fname)
    
    #create the HDF5 file
    try:
        f1 = h5py.File(fname.rsplit('.mrc',1)[0] + '.emd','w') #w- will error if the file exists
    except:
        print("Problem opening file. Maybe it already exists?")
        f1.close()
        del tomo
        return 0

    #Create the axis vectors in nanometers. Standard MRC pixel size is in Angstroms
    xFull = np.linspace(0,tomo['voxelSize'][0]*tomo['stack'].shape[0]-1,tomo['stack'].shape[0]) 
    yFull = np.linspace(0,tomo['voxelSize'][1]*tomo['stack'].shape[1]-1,tomo['stack'].shape[1])
    zFull = np.linspace(0,tomo['voxelSize'][2]*tomo['stack'].shape[2]-1,tomo['stack'].shape[2])
    
    #Root data group
    dataTop = f1.create_group('data')

    #Create tilt series group
    tiltseriesGroup = dataTop.create_group('stack')
    tiltseriesGroup.attrs['emd_group_type'] = np.int8(1)
    
    #Save the data to the EMD file and reshape it to a C-style array
    try:
        tiltDset = tiltseriesGroup.create_dataset('data',data=tomo['stack'],compression='gzip',shuffle=True)
    except MemoryError:
        print("Not enough memory to write out data to EMD file")
        del tomo
        f1.close()
        return 0
        
    dim1 = tiltseriesGroup.create_dataset('dim1',data=xFull)
    dim1.attrs['name'] = np.string_('x')
    dim1.attrs['units'] = np.string_('')
    dim2 = tiltseriesGroup.create_dataset('dim2',data=yFull)
    dim2.attrs['name'] = np.string_('y')
    dim2.attrs['units'] = np.string_('')
    dim3 = tiltseriesGroup.create_dataset('dim3',data=zFull)
    dim3.attrs['name'] = np.string_('z')
    dim3.attrs['units'] = np.string_('')
    
    #Create the other groups
    scopeGroup = f1.create_group('Microscope')
    scopeGroup.attrs['voxel sizes'] = tomo['voxelSize']
    userGroup = f1.create_group('User')
    commentGroup = f1.create_group('Comments')
    
    #Possible way using keyword arguments to populate these fields
    #def greet_me(**kwargs):
    #if kwargs is not None:
    #    for key, value in kwargs.iteritems():
    #        print("%s == %s" %(key,value))
    
    f1.close()
    
    return 1
    
def mrcWriter(filename,stack,pixelSize,forceWrite=False):
    '''Write out a MRC type file according to the specification at http://bio3d.colorado.edu/imod/doc/mrc_format.txt
    
    Parameters:
        filename (str): The name of the MRC file.
        stack (ndarray): The array data to write to disk.
        pixelSize (tuple): The size of the pixel along each direction (in Angstroms) as a 3 element vector (sizeZ,sizeY,sizeX).
    Returns:
        (int): 1 if successful and 0 if unsuccessful
    '''
    
    fid = open(filename,'wb')
    
    if len(stack.shape) > 3:
        print("Too many dimensions")
        return 0;
    
    if not stack.flags['C_CONTIGUOUS']:
        print("Error: Array must be C-style ordering: [numImages,Y,X]. Use numpy.tranpspose and np.ascontiguousarray to change data ordering in memory")
        print('Exiting')
        return 0;
    
    #initialize the header with 256 zeros with size 4 bytes
    header = np.zeros(256,dtype=np.int32)
    fid.write(header)
    fid.seek(0,0) #return to the beginning of the file
    
    #Initialize the int32 part of the header
    header1 = np.zeros(10,dtype=np.int32)
    
    #Write the number of columns, rows and sections (images)
    #header1[0:3] = np.int32(dims) #stack size in pixels
    header1[0] = np.int32(stack.shape[2]) #num columns, the last index in C-style ordering
    header1[1] = np.int32(stack.shape[1]) #num rows
    header1[2] = np.int32(stack.shape[0]) #num sections (images)
    
    if stack.dtype == np.float32:
        header1[3] = np.int32(2)
    elif stack.dtype == np.uint16:
        header1[3] = np.int32(6)
    elif stack.dtype == np.int16:
        header1[3] = np.int32(1)
    elif stack.dtype == np.int8:
        header1[3] = np.int32(0)
    else:
        print("Data type " + str(stack.dtype) + " is unsupported. Only int8, int16, uint16, and float32 are supported")
        return 0;
    
    #Starting point of sub image (not used in IMOD) 
    header1[4:7] = np.zeros(3,dtype=np.int32)
    
    #Grid size in X,Y,Z
    #header1[7:10] = np.int32(dims); #stack size in pixels
    header1[7] = np.int32(stack.shape[2]) #mx
    header1[8] = np.int32(stack.shape[1]) #my
    header1[9] = np.int32(stack.shape[0]) #mz
    
    #Write out the first part of the header information
    fid.write(header1)
    
    #Cell dimensions (in Angstroms)
    #pixel spacing = xlen/mx, ylen/my, zlen/mz
    fid.write(np.float32(pixelSize[2]*stack.shape[2])) #xlen
    fid.write(np.float32(pixelSize[1]*stack.shape[1])) #ylen
    fid.write(np.float32(pixelSize[0]*stack.shape[0])) #zlen
    
    #Cell angles (in degrees)
    fid.write(np.float32([90.0,90.0,90.0]))
    
    #Description of array directions with respect to: Columns, Rows, Images
    fid.write(np.int32([1,2,3]))
    
    #Minimum and maximum density
    fid.write(np.float32(np.min(stack)))
    fid.write(np.float32(np.max(stack)))
    fid.write(np.float32(np.mean(stack)))
    
    #Needed to indicate that the data is little endian for NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above
    fid.seek(212,0)
    fid.write(np.int8([68,65,0,0])) #use [17,17,0,0] for big endian
    
    #Write out the data
    fid.seek(1024)
    if forceWrite:    
        fid.write(np.ascontiguousarray(stack)) #Change to C ordering array for writing to disk
    else:
        fid.write(stack) #msut be C-contiguous

    #Close the file
    fid.close()
    return 1

def writeHeader(filename,shape,dtype,pixelSize):
    '''Write out a MRC type file header according to the specification at http://bio3d.colorado.edu/imod/doc/mrc_format.txt.
    This is useful for initializing an MRC file and then writing to it manually or see appendData() function below.
    
    Parameters:
        filename (str): The name of the EMD file 
        shape (tuple): The shape of the data to write
        pixelSize (tuple): The size of the pixel along each direction (in Angstroms) as a 3 element vector (sizeX,sizeY,sizeZ). sizeZ could be the angular step for a tilt series
    Returns:
        (int): 1 if successful.
    '''
    
    with open(filename,'wb') as fid:
        if len(shape) > 3:
            print("Too many dimensions")
            return 0
        
        #if not stack.flags['C_CONTIGUOUS']:
        #    print("Error: Array must be C-style ordering: [numImages,Y,X]. Use np.ascontiguousarray(numpy.transpose(stack,[])) to change data ordering in memory\nExiting")
        #    return 0
        
        #Initialize the header with 256 zeros with size 4 bytes
        header = np.zeros(256,dtype=np.int32)
        fid.write(header)
        fid.seek(0,0) #return to the beginning of the file
        
        #Initialize the int32 part of the header with zeros
        header1 = np.zeros(10,dtype=np.int32)
        
        #Write the number of columns, rows and sections (images)
        header1[0] = np.int32(shape[2]) #num columns, the last index in C-style ordering
        header1[1] = np.int32(shape[1]) #num rows
        header1[2] = np.int32(shape[0]) #num sections (images)
        
        if dtype == np.float32:
            header1[3] = np.int32(2)
        elif dtype == np.uint16:
            header1[3] = np.int32(6)
        elif dtype == np.int16:
            header1[3] = np.int32(1)
        elif dtype == np.int8:
            header1[3] = np.int32(0)
        else:
            print("Data type " + str(dtype) + " is unsupported. Only int8, int16, uint16, and float32 are supported")
            return 0;
        
        #Starting point of sub image (not used in IMOD) 
        header1[4:7] = np.zeros(3,dtype=np.int32)
        
        #Grid size in X,Y,Z
        header1[7] = np.int32(shape[2]) #mx
        header1[8] = np.int32(shape[1]) #my
        header1[9] = np.int32(shape[0]) #mz
        
        #Write out the first part of the header information
        fid.write(header1)
        
        #Cell dimensions (in Angstroms)
        #pixel spacing = xlen/mx, ylen/my, zlen/mz
        fid.write(np.float32(pixelSize[2]*shape[2])) #xlen
        fid.write(np.float32(pixelSize[1]*shape[1])) #ylen
        fid.write(np.float32(pixelSize[0]*shape[0])) #zlen
        
        #Cell angles (in degrees)
        fid.write(np.float32([90.0,90.0,90.0]))
        
        #Description of array directions with respect to: Columns, Rows, Images
        fid.write(np.int32([1,2,3]))
        
        #Minimum and maximum density
        #fid.write(np.float32(np.min(stack)))
        #fid.write(np.float32(np.max(stack)))
        #fid.write(np.float32(np.mean(stack)))
        fid.write(np.zeros(3,dtype=np.float32))
        
        #Needed to indicate that the data is little endian for NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above
        fid.seek(212,0)
        fid.write(np.int8([68,65,0,0])) #use [17,17,0,0] for big endian
        
    return 1

def appendData(filename,stack):
    '''Append a binary set of data to the end of a MRC file. This should only be used in conjunction with
    writeHeader() above. 
    
    Parameters:
        filename (str):    Name of the MRC file with pre-initiated header and some data already written.
        stack (ndarray):   Data to append to the file.
    
    '''
    with open(filename,'ab') as fid:
        #Write out the data
        fid.seek(0,2) #seek to the end of the file
        fid.write(stack) #Change to C ordering array for writing to disk
        
def emd2mrc(filename,dsetPath):
    '''Convert EMD data set into MRC data set. The final data type is float32 for convenience.
    
    Parameters:
        filename (str): The name of the EMD file
        dsetPath (str): The HDF5 path to the top group holding the data. ex. '/data/raw/'
    '''
    
    with h5py.File(filename,'r') as f1:
        #Get the pixel sizes and convert to Ang
        dimsPath = dsetPath + '/dim'
        print('Warning: Assuming EMD dim vectors are in nanometer')
        print('dim vector names/units are:')
        for ii in range(1,4):
            print('name, units = {}, {}'.format(f1[dimsPath+str(ii)].attrs['name'],f1[dimsPath+str(ii)].attrs['units']))
        pixelSizeX = (f1[dsetPath + '/dim2'][1] - f1[dsetPath + '/dim2'][0])*10 #change nanometers to Ang
        pixelSizeY = (f1[dsetPath + '/dim3'][1] - f1[dsetPath + '/dim3'][0])*10 #change nanometers to Ang
                
        filenameOut = filename.split('.emd')[0] + '.mrc' #use the first part of the file as the prefix removing the .emd on the end
        
        print('Warning: Converting to float32 before writing to disk')
        mrc.mrcWriter(filenameOut,np.float32(f1[dsetPath+'/data']),(1,pixelSizeY,pixelSizeX)) #the extra slash is not a problem. // is the same as / in a HDF5 data set path
        
        print('Finished writing to: {}'.format(filenameOut))