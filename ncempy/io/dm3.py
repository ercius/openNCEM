"""
A module to load data and meta data from DM3 files into python

"""
import numpy as np

class fileDM3:
    def __init__(self, filename, verbose = False):
        '''Init opening the file and reading in the header.
        '''
        # necessary declarations, if something fails
        self.fid = None
        self.fidOut = None
        
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string')
        
        #Add a top level variable to indicate verbosee output for debugging
        self.v = verbose
        
        # try opening the file
        try:
            self.fid = open(filename, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(filename))
            raise
        except :
            raise
        
        if not self.validDM3():
            print('Not a valid DM3 file: "{}"'.format(filename))
            raise IOError('Improperly formated file')
        
        #Lists that will contain information about binary data arrays
        self.xSize = []
        self.ySize = []
        self.zSize = []
        self.dataType = []
        self.dataSize = []
        self.dataOffset = []
        
        #The number of objects found in the DM3 file
        self.numObjects = 0
        
        self.curGroupLevel = 0 #track how deep we currently are in a group
        self.maxDepth = 64 #maximum number of group levels allowed
        self.curGroupAtLevelX = np.zeros((self.maxDepth,),dtype=np.int8) #track group at current level
        self.curGroupNameAtLevelX = '' #set the name of the root group
        
        self.curTagAtLevelX = np.zeros((self.maxDepth,),dtype=np.int8) #track tag number at the current level
        self.curTagName = '' #string of the current tag
        
        #lists that will contain scale information (pixel size)
        self.scale = [] 
        self.scaleUnit = []
        self.origin = []
        
        #Temporary variables to keep in case a tag entry shows useful information in an array
        self.scale_temp = 0
        self.origin_temp = 0
        
        self.outputDic = {}
        self.allTags = {}
        
        #Read the DM3 header as a set of DM tags
        self.readTagGroup()
    
    def __del__(self):
        #close the file
        if(self.fid):
            print('Closing input file')
            self.fid.close()
        if(self.fidOut):
            print('Closing output file')
            self.fidOut.close()
    
    def validDM3(self):
        '''Test whether a file is a valid DM3 file and written in Little Endian format
        '''
        output = True #output will stay == 1 if the file is a true DM3 file

        head = np.fromfile(self.fid,dtype=np.dtype('>u4'),count=3)
        
        if head[0] != 3:
            print('File is not a dm3. DM file type number is {}'.format(head[0]))
            output = False
        
        #Useful to test against current file size
        self.fileSize = head[1]
            
        if head[2] != 1:
            print('File is not written Little Endian (PC) format and can not be read by this program.')
            output = False
        
        return output
    
    def readTagGroup(self):
        
        self.curGroupLevel += 1
        self.curGroupAtLevelX[self.curGroupLevel] = self.curGroupAtLevelX[self.curGroupLevel] + 1
        self.curTagAtLevelX[self.curGroupLevel] = 0
        np.fromfile(self.fid,dtype='<i1',count=2) #is open and is sorted?
        nTags = np.fromfile(self.fid,dtype='>u4',count=1)[0] #needs to be read as Big Endian (.byteswap() could also work)
        
        #Iterate of the number of tag entries
        oldTotalTag = self.curGroupNameAtLevelX
        for ii in range(0,nTags):
            self.readTagEntry()
        
        #Go back up a level after reading all entries
        self.curGroupLevel -= 1
        self.curGroupNameAtLevelX = oldTotalTag
        
    def readTagEntry(self):
        dataType = np.fromfile(self.fid,dtype=np.dtype('>u1'),count=1)[0]
        
        if self.v:
            print('readTagEntry: dataType = {}'.format(dataType))
        
        #Record tag at this level
        self.curTagAtLevelX[self.curGroupLevel] += 1
        
        #get the tag
        lenTagLabel = np.fromfile(self.fid,dtype=np.dtype('>u2'),count=1)[0]
        if lenTagLabel > 0:
            tagLabelBinary = np.fromfile(self.fid,dtype='<u1',count=lenTagLabel) #read as binary
            #tagLabel = ''.join([chr(item) for item in tagLabelBinary]) #convert to python string
            tagLabel = self.bin2str(tagLabelBinary)
        else:
            tagLabel = str(self.curTagAtLevelX[self.curGroupLevel]) #unlabeled tag.
        
        #Save the current group name in case this is needed
        oldGroupName = self.curGroupNameAtLevelX
        
        if dataType == 21:
            #This tag entry contains data
            self.curTagName = tagLabel #save its name
            self.readTagType()
        else:
            #This is a nested tag group
            self.curGroupNameAtLevelX += '.' + tagLabel #add to group names
            self.readTagGroup()
        self.curGroupNameAtLevelX = oldGroupName
    
    def readTagType(self):
        delim = np.fromfile(self.fid,dtype='<i1',count=4)
        assert((delim == 37).all()) #delim has to be [37,37,37,37] which is %%%% in ASCII.
        if self.v:
            print('readTagType: should be %%%% = {}'.format(self.bin2str(delim)))
        nInTag = np.fromfile(self.fid,dtype=np.dtype('>u4'),count=1)[0] #nInTag: unnecessary redundant info

        #readAnyData() #calling this separately is unnecessary
     
     #def readAnyData() #remove separate readAnyData
        
        #Determine the type of the data in the tag
        #specifies data type: int8, uint16, float32, etc.
        encodedType = np.fromfile(self.fid,dtype=np.dtype('>u4'),count=1)[0] #big endian
        
        etSize = self.encodedTypeSize(encodedType)
        
        if etSize > 0:
            #regular data. Read it and store it with the tag name
            if self.v:
                print('regular')
            self.storeTag(self.curTagName, self.readNativeData(encodedType))
        elif encodedType == 18: #string
            if self.v:
                print('string')
            stringSize = np.fromfile(self.fid,dtype='>u4',count=1)[0]
            #strtemp = '' #in case stringSize == 0
            strTempBin = np.fromfile(self.fid,dtype='<u1',count=stringSize) #read as uint8 little endian
            strTemp = self.bin2str(strTempBin)
            self.storeTag(self.curTagName,strTemp)
        elif encodedType == 15: #struct
            #This does not work for field names that are non-zero. This is uncommon
            if self.v:
                print('struct')
            structTypes = self.readStructTypes()
            structs = self.readStructData(structTypes)
            self.storeTag(self.curTagName,structs)
        elif encodedType == 20: #array
            #The array data is not read. It will be read later if needed
            if self.v:
                print('array')
            arrayTypes = self.readArrayTypes() #could be recursive if array contains array(s)
            arrInfo = self.readArrayData(arrayTypes) #only info of the array is read. It is read later if needed
            self.storeTag(self.curTagName,arrInfo)
    
    def bin2str(self,bin):
        '''Utility function to convert an numpy array of binary values to a python string
        '''
        if self.v:
            print('bin2str: input binary numbers = {}'.format(bin))
        return ''.join([chr(item) for item in bin])
        
    def encodedTypeSize(self, encodedType):
        '''Return the number of bytes in a data type for the encodings used by DM
        Constants for the different encoded data types used in DM3 files
         	VAL_SHORT   = 2;
         	VAL_LONG    = 3;
         	VAL_USHORT  = 4;
         	VAL_ULONG   = 5;
         	VAL_FLOAT   = 6;
         	VAL_DOUBLE  = 7;
         	VAL_BOOLEAN = 8;
         	VAL_CHAR    = 9;
         	VAL_OCTET   = 10;
           -1 will signal an unlisted type
        '''
        if encodedType == 0:
            return 0
        elif (encodedType == 8) | (encodedType == 9) | (encodedType == 10):
            return 1 #1 byte
        elif (encodedType == 2) | (encodedType == 4):
            return 2 #2 bytes
        elif (encodedType == 3) | (encodedType == 5) | (encodedType == 6):
            return 4 #4 bytes
        elif (encodedType == 7):
            return 8 #8 bytes
        else:
            return -1
    
    def readStructTypes(self):
        '''Analyze the types of data in a struct
        '''
        structNameLength = np.fromfile(self.fid,count=1,dtype='>u4')[0] #this is not needed
        nFields = np.fromfile(self.fid,count=1,dtype='>u4')[0]
        if self.v:
            print('readStructTypes: nFields = {}'.format(nFields))
        
        if(nFields > 100):
            raise RuntimeError('Too many fields in a struct.')
        
        fieldTypes = np.zeros(nFields)
        for ii in range(0,nFields):
            aa = np.fromfile(self.fid,count=2,dtype=np.dtype('>u4')) #nameLength, fieldType
            fieldTypes[ii] = aa[1]
        return fieldTypes
    
    def readStructData(self,structTypes):
        '''Read the data in a struct
        '''
        struct = np.zeros(structTypes.shape[0])
        for ii, encodedType in enumerate(structTypes):
            etSize = self.encodedTypeSize(encodedType) #the size of the data type
            struct[ii] = self.readNativeData(encodedType) #read this type of data
        return struct
    
    def readNativeData(self,encodedType):
        #reads ordinary data types
        # 	VAL_SHORT   = 2;
        # 	VAL_LONG    = 3;
        # 	VAL_USHORT  = 4;
        # 	VAL_ULONG   = 5;
        # 	VAL_FLOAT   = 6;
        # 	VAL_DOUBLE  = 7;
        # 	VAL_BOOLEAN = 8;
        # 	VAL_CHAR    = 9;
        # 	VAL_OCTET   = 10;   
        #   VAL_UINT64 = 11;
        if encodedType == 2:
            val = np.fromfile(self.fid,count=1,dtype='<i2')[0]
        elif encodedType == 3:
            val = np.fromfile(self.fid,count=1,dtype='<i4')[0]
        elif encodedType == 4:
            val = np.fromfile(self.fid,count=1,dtype='<u2')[0]
        elif encodedType == 5:
            val = np.fromfile(self.fid,count=1,dtype='<u4')[0]
        elif encodedType == 6:
            val = np.fromfile(self.fid,count=1,dtype='<f4')[0]
        elif encodedType == 7:
            val = np.fromfile(self.fid,count=1,dtype='<f8')[0]
        elif encodedType == 8: #matlab uchar
            print('readNativeData untested type: {}'.format(encodedType))
            val = np.fromfile(self.fid,count=1,dtype='<u1')[0] #return character or number?
        elif encodedType == 9: #matlab *char
            print('readNativeData untested type: {}'.format(encodedType))
            val = np.fromfile(self.fid,count=1,dtype='<i1')[0] #return character or number?
        elif encodedType == 10: #matlab *char
            print('readNativeData untested type: {}'.format(encodedType))
            val = np.fromfile(self.fid,count=1,dtype='<i1')[0]
        elif encodedType == 11:
            val = np.fromfile(self.fid,count=1,dtype='<i8')[0]
        else:
            print('readNativeData unknown data type: {}'.format(encodedType))
            raise
        
        if self.v:
            print('readNativeData: encodedType == {} and val = {}'.format(encodedType, val))
        
        return val
    def readArrayTypes(self):
        '''Analyze te types of data in an array
        '''
        arrayType = np.fromfile(self.fid,count=1,dtype='>u4')[0]
        
        itemTypes = []
        
        if arrayType == 15: 
            #nested Struct
            itemTypes = self.readStructTypes()
        elif arrayType == 20:
            #Nested array
            itemTypes = readArrayTypes()
        else:
            itemTypes.append(arrayType)
        if self.v:
            print('readArrayTypes: itemTypes = {}'.format(itemTypes))
        return itemTypes
    
    def readArrayData(self,arrayTypes):
        '''Read information in an array based on the types provided. Binary data is not read at this point.
        '''
        arraySize = np.fromfile(self.fid,count=1,dtype='>u4')[0]
        
        itemSize = 0
        encodedType = 0
        
        if self.v:
            print('readArrayData: arrayTypes = {}'.format(arrayTypes))
        
        for encodedType in arrayTypes:
            print('encodedType = {}'.format(encodedType))
            etSize = self.encodedTypeSize(encodedType)
            itemSize += etSize
            
        bufSize = arraySize * itemSize
        
        print('arraySize, itemSize = {}, {}'.format(arraySize, itemSize))
        
        if encodedType == 4:
            #String data
            self.storeTag(self.curTagName + '.arraySize',bufSize)
            self.storeTag(self.curTagName + '.arrayOffset', self.fid.tell())
            stringData = self.bin2str(np.fromfile(self.fid,count=bufSize,dtype='<u1'))
            #print('readArrayData: stringData =  {}'.format(stringData))
            arrOut = stringData.replace('\x00','') #remove all spaces from the string data
            #print('readArrayData: arrOut =  {}'.format(arrOut))
            
            #Catch useful tags for images and spectra (nm, eV, etc.)
            fullTagName = self.curGroupNameAtLevelX + '.' + self.curTagName
            if((fullTagName.find('Dimension') > -1) & (fullTagName.find('Units') > -1) & (self.numObjects > 0)):
                #tt = len(self.scale) #the number of dimension and units already saved
                self.scale.append(self.scale_temp)
                self.scaleUnit.append(arrOut)
                self.origin.append(self.origin_temp)
        else:
            #This is a binary array. Save its location to read later if needed
            self.storeTag(self.curTagName + '.arraySize', bufSize)
            self.storeTag(self.curTagName + '.arrayOffset', self.fid.tell())
            self.storeTag(self.curTagName + '.arrayType', encodedType)
            self.fid.seek(bufSize,1) #advance the pointer by bufsize from current position
            arrOut = 'Array data unread. Encoded type = {}'.format(encodedType)
            
        return arrOut
    
    def storeTag(self,curTagName,curTagValue):
        '''Builds the full tag name and key/value pair as text. Also calls another
        function to catch useful tags and values. Also saves all tags in a dictionary.
        '''
        #Build the full tag name (key) and add the tag value
        if self.v:
            print('storeTag: curTagName, curTagValue = {}, {}'.format(curTagName,curTagValue))
        totalTag = self.curGroupNameAtLevelX + '.' + '{}'.format(curTagName) + '= {}'.format(curTagValue)
        
        self.catchUsefulTags(totalTag,curTagName,curTagValue)
        
        self.allTags[totalTag] = curTagValue #this needs to be done better. 
        
        return(totalTag)
    
    def catchUsefulTags(self,totalTag,curTagName,curTagValue):
        '''Find interesting keys and keep their values for later. This is separate from storeTag
        so that it is easy to find and modify.
        '''
        if curTagName.find('Data.arraySize')>-1:
            self.numObjects += 1 #add this as an interesting object
            self.dataSize.append(curTagValue)
        elif curTagName.find('Data.arrayOffset') >-1:
            self.dataOffset.append(curTagValue)
        elif curTagName.find('DataType')>-1:
            self.dataType.append(curTagValue)
        elif totalTag.find('Dimensions.1')>-1:
            self.xSize.append(curTagValue)
            self.ySize.append(1) 
            self.zSize.append(1)
        elif totalTag.find('Dimensions.2')>-1:
            self.ySize[-1] = curTagValue #OR self.ysize[self.numObjects] = self.curTagValue
        elif totalTag.find('Dimensions.3')>-1:
            self.zSize[-1] = curTagValue
        elif (totalTag.find('Dimension.')>-1) & (totalTag.find('.Scale')>-1):
            self.scale_temp = curTagValue
        elif (totalTag.find('Dimension.')>-1) & (totalTag.find('.Origin')>-1):
            self.origin_temp = curTagValue
        else:
            pass
    
    def writeTags(self):
        fnameOutPrefix = self.filename.split('.dm3')[0]
        try:
            #open a text file to write out the tags
            with open(fnameOutPrefix+'_tags.txt','w') as fidOut:
                for nn in self.allTags:
                    fidOut.write(nn + ' = ' + str(self.allTags[nn]))
            fidOut.close() #this might not be necessary
        except NameError:
            print("Issue opening tags output file.")
            raise
        except:
            raise
    
    def checkIndex(self, i):
        '''Check index i for sanity, otherwise raise Exception.
        
        Parameters:
            i (int):    Index.
            
        '''
        
        # check type
        if not isinstance(i, int):
            raise TypeError('index supposed to be integer')

        # check whether in range
        if i < 0 or i > self.numObjects:
            raise IndexError('Index out of range, trying to access element {} of {} valid elements'.format(i+1, self.head['ValidNumberElements']))
            
        return        
    
    def getNPDataType(self, dd):
        '''Convert the DM data type value into a numpy datatype
        '''
        if dd == 6:
            return np.uint8
        elif dd == 10:
            return np.uint16
        elif dd == 11:
            return np.uint32
        elif dd == 9:
            return np.int8
        elif dd == 1:
            return np.int16
        elif dd == 7:
            return np.int32
        elif dd == 2:
            return np.float32
        elif dd == 12:
            return np.float64
        #elif dd == 14: #this is supposed to be bit1 in matlab, but Im not sure what that translates to in numpy
        #    return np.uint8 #bit1 ??
        elif dd == 3:
            return np.complex64
        elif dd == 13:
            return np.complex128
        else:
            print('dm3reader: Unsupported data type: DM dataType == {}'.format(dd))
    
    def getDataset(self, index):
        '''Retrieve a dataseet from the DM3 file.
        '''
        
        #The first dataset is always a thumbnail. Skip it
        ii = index + 1
        
        #Check that the dataset exists.
        try:
            self.checkIndex(ii)
        except:
            raise
        
        self.fid.seek(self.dataOffset[ii],0) #Seek to start of dataset from beginning of the file
        
        outputDict = {}
        
        #Parse the dataset to see what type it is (image, image series, spectra, etc.)
        if self.xSize[ii] > 0:
            outputDict['scaleUnit'] = self.scaleUnit[ii]
            pixelCount = self.xSize[ii]*self.ySize[ii]*self.zSize[ii]
            #if self.dataType == 23: #RGB image(s)
            #    temp = np.fromfile(self.fid,count=pixelCount,dtype=np.uint8).reshape(self.ysize[ii],self.xsize[ii])
            if self.zSize[ii] == 1:
                outputDict['image'] = np.fromfile(self.fid,count=pixelCount,dtype=self.getNPDataType(self.dataType[ii])).reshape((self.ySize[ii],self.xSize[ii]))
            else: #3D array
                outputDict['cube'] = np.fromfile(self.fid,count=pixelCount,dtype=self.getNPDataType(self.dataType[ii])).reshape((self.zSize[ii],self.ySize[ii],self.xSize[ii]))
                #outputDict['cube'] = np.fromfile(self.fid,count=pixelCount,dtype=np.int16).reshape((self.zSize[ii],self.ySize[ii],self.xSize[ii]))
            
        
        return outputDict