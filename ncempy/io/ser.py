'''
This module provides an interface to the SER file format written by TIA.

Following the information provided by Dr Chris Boothroyd (http://www.er-c.org/cbb/info/TIAformat/).
'''
import numpy as np
import h5py
import os
import re
import xml.etree.ElementTree as ET
import datetime
import ncempy.io.emd

class NotSERError(Exception):
    '''Exception if a file is not in SER file format.
    
    '''
    pass


class fileSER:
    '''
    Class to represent SER files (read only).
    
    Parameters:
        filename (str):    Name of the SER file.
        emifile (str):    Name of an optional _emi file to read metadata.
        verbose (bool):    True to get extensive output while reading the file.
        
    '''

    _dictByteOrder = {0x4949 : 'little endian'}
    '''(dict):    Information on byte order.'''
    
    _dictSeriesVersion = {0x0210 : '< TIA 4.7.3', 0x0220 : '>= TIA 4.7.3'}
    '''(dict):    Information on fileformat version.'''
    
    _dictDataTypeID = {0x4120:'1D datasets', 0x4122:'2D images'}
    '''(dict):    Information on data type.'''
    
    _dictTagTypeID = {0x4152:'time only',0x4142:'time and 2D position'}
    '''(dict):    Information on tag type.'''
    
    _dictDataType = {1:'<u1', 2:'<u2', 3:'<u4', 4:'<i1', 5:'<i2', 6:'<i4', 7:'<f4', 8:'<f8', 9:'<c8', 10:'<c16'}
    '''(dict):    Information on data format.'''


    def __init__(self, filename, emifile=None, verbose=False):
        '''Init opening the file and reading in the header.
        
        '''
        # necessary declarations, if something fails
        self._file_hdl = None
        self._emi = None

        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string')

        # try opening the file
        try:
            self._file_hdl = open(filename, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(filename))
            raise
        except :
            raise

        # read header
        self.head = self.readHeader(verbose)
        
        # read _emi, if provided
        if emifile:
            self._emi = self.read_emi(emifile)


    def __del__(self):
        '''Closing the file stream on del.
        
        '''
        
        # close the file
        if(self._file_hdl):
            self._file_hdl.close()


    def readHeader(self, verbose=False):
        '''Read and return the SER files header.
        
        Parameters:
            verbose (bool):    True to get extensive output while reading the file.

        Returns:
            (dict):    The header of the SER file as dict.
            
        '''
        
        # prepare empty dict to be populated while reading
        head = {}

        # go back to beginning of file
        self._file_hdl.seek(0,0)

        # read 3 int16
        data = np.fromfile(self._file_hdl, dtype='<i2', count=3)

        # ByteOrder (only little Endian expected)
        if not data[0] in self._dictByteOrder:
            raise RuntimeError('Only little Endian implemented for SER files')
        head['ByteOrder'] = data[0]
        if verbose:
            print('ByteOrder:\t"{:#06x}",\t{}'.format(data[0], self._dictByteOrder[data[0]]))

        # SeriesID, check whether TIA Series Data File   
        if not data[1] == 0x0197:
            raise NotSERError('This is not a TIA Series Data File (SER)')
        head['SeriesID'] = data[1]
        if verbose:
            print('SeriesID:\t"{:#06x},\tTIA Series Data File'.format(data[1]))

        # SeriesVersion
        if not data[2] in self._dictSeriesVersion:
            raise RuntimeError('Unknown TIA version: "{:#06x}"'.format(data[2]))
        head['SeriesVersion'] = data[2]
        if verbose:
            print('SeriesVersion:\t"{:#06x}",\t{}'.format(data[2], self._dictSeriesVersion[data[2]]))
        # version dependend fileformat for below
        if head['SeriesVersion']==0x0210:
            offset_dtype = '<i4'
        else: #head['SeriesVersion']==0x220:
            offset_dtype = '<i8'
        
        # read 4 int32
        data = np.fromfile(self._file_hdl, dtype='<i4', count=4)

        # DataTypeID
        if not data[0] in self._dictDataTypeID:
            raise RuntimeError('Unknown DataTypeID: "{:#06x}"'.format(data[0]))
        head['DataTypeID'] = data[0]
        if verbose:
            print('DataTypeID:\t"{:#06x}",\t{}'.format(data[0], self._dictDataTypeID[data[0]]))

        # TagTypeID
        if not data[1] in self._dictTagTypeID:
            raise RuntimeError('Unknown TagTypeID: "{:#06x}"'.format(data[1]))
        head['TagTypeID'] = data[1]
        if verbose:
            print('TagTypeID:\t"{:#06x}",\t{}'.format(data[1], self._dictTagTypeID[data[1]]))

        # TotalNumberElements
        if not data[2] >= 0:
            raise RuntimeError('Negative total number of elements: {}'.format(data[2]))
        head['TotalNumberElements'] = data[2]
        if verbose:
            print('TotalNumberElements:\t{}'.format(data[2]))

        # ValidNumberElements
        if not data[3] >= 0:
            raise RuntimeError('Negative valid number of elements: {}'.format(data[3]))
        head['ValidNumberElements'] = data[3]
        if verbose:
            print('ValidNumberElements:\t{}'.format(data[3]))
        
        # OffsetArrayOffset, sensitive to SeriesVersion
        data = np.fromfile(self._file_hdl, dtype=offset_dtype, count=1)
        head['OffsetArrayOffset'] = data[0]
        if verbose:
            print('OffsetArrayOffset:\t{}'.format(data[0]))

        # NumberDimensions
        data = np.fromfile(self._file_hdl, dtype='<i4', count=1)
        if not data[0] >= 0:
            raise RuntimeError('Negative number of dimensions')
        head['NumberDimensions']=data[0]
        if verbose:
            print('NumberDimensions:\t{}'.format(data[0]))


        # Dimensions array
        dimensions = []
        for i in range(head['NumberDimensions']):
            if verbose:
                print('reading Dimension {}'.format(i))
            this_dim = {}
        
            # DimensionSize
            data = np.fromfile(self._file_hdl, dtype='<i4', count=1)
            this_dim['DimensionSize'] = data[0]
            if verbose:
                print('DimensionSize:\t{}'.format(data[0]))
            
            data = np.fromfile(self._file_hdl, dtype='<f8', count=2)
            
            # CalibrationOffset
            this_dim['CalibrationOffset'] = data[0]
            if verbose:
                print('CalibrationOffset:\t{}'.format(data[0]))
            
            # CalibrationDelta
            this_dim['CalibrationDelta'] = data[1]
            if verbose:
                print('CalibrationDelta:\t{}'.format(data[1]))
            
            data = np.fromfile(self._file_hdl, dtype='<i4', count=2)
            
            # CalibrationElement
            this_dim['CalibrationElement'] = data[0]
            if verbose:
                print('CalibrationElement:\t{}'.format(data[0]))
            
            # DescriptionLength
            n = data[1]
            
            # Description
            data = np.fromfile(self._file_hdl, dtype='<i1', count=n)
            data = ''.join(map(chr, data))
            this_dim['Description'] = data
            if verbose:
                print('Description:\t{}'.format(data))
            
            # UnitsLength
            data = np.fromfile(self._file_hdl, dtype='<i4', count=1)
            n = data[0]
            
            # Units
            data = np.fromfile(self._file_hdl, dtype='<i1', count=n)
            data = ''.join(map(chr, data))
            this_dim['Units'] = data
            if verbose:
                print('Units:\t{}'.format(data))


            dimensions.append(this_dim)
        
        # save dimensions array as tuple of dicts in head dict
        head['Dimensions'] = tuple(dimensions)
        
        
        # Offset array
        self._file_hdl.seek(head['OffsetArrayOffset'],0)
        
        # DataOffsetArray
        data = np.fromfile(self._file_hdl, dtype=offset_dtype, count=head['ValidNumberElements'])
        head['DataOffsetArray'] = data.tolist()
        if verbose:
            print('reading in DataOffsetArray')
        
        # TagOffsetArray
        data = np.fromfile(self._file_hdl, dtype=offset_dtype, count=head['ValidNumberElements'])
        head['TagOffsetArray'] = data.tolist()
        if verbose:
            print('reading in TagOffsetArray')     

        return head


    def _checkIndex(self, i):
        '''Check index i for sanity, otherwise raise Exception.
        
        Parameters:
            i (int):    Index.
            
        '''
        
        # check type
        if not isinstance(i, int):
            raise TypeError('index supposed to be integer')

        # check whether in range
        if i < 0 or i>= self.head['ValidNumberElements']:
            raise IndexError('Index out of range, trying to access element {} of {} valid elements'.format(i+1, self.head['ValidNumberElements']))
            
        return
        

    def getDataset(self, index, verbose=False):
        '''Retrieve dataset from data file.

        Parameters:
            index (int):    Index of dataset.
            verbose (bool):    True to get extensive output while reading the file.
        
        Returns:
            (tuple):    Tuple containing:
            
                np.ndarray:    Dataset as array.
                
                dict:    Metadata as dict.
                
        '''

        # check index, will raise Exceptions if not
        try:
            self._checkIndex(index)
        except:
            raise

        if verbose:
            print('Getting dataset {} of {}.'.format(index, self.head['ValidNumberElements']))
            
        # go to dataset in file
        self._file_hdl.seek(self.head['DataOffsetArray'][index],0)
        
        # read meta
        meta = {}
        
        # number of calibrations depends on DataTypeID
        if self.head['DataTypeID'] == 0x4120:
            n = 1
        elif self.head['DataTypeID'] == 0x4122:
            n = 2
        else:
            raise RuntimeError('Unknown DataTypeID')
       
        # read in the calibrations    
        cals = []
        for i in range(n):
            if verbose:
                print('Reading calibration {}'.format(i))
                
            this_cal = {}
        
            data = np.fromfile(self._file_hdl, dtype='<f8', count=2)
            
            # CalibrationOffset
            this_cal['CalibrationOffset'] = data[0]
            if verbose:
                print('CalibrationOffset:\t{}'.format(data[0]))
            
            # CalibrationDelta
            this_cal['CalibrationDelta'] = data[1]
            if verbose:
                print('CalibrationDelta:\t{}'.format(data[1]))
            
            data = np.fromfile(self._file_hdl, dtype='<i4', count=1)
            
            # CalibrationElement
            this_cal['CalibrationElement'] = data[0]
            if verbose:
                print('CalibrationElement:\t{}'.format(data[0]))
            
            cals.append(this_cal)
            
        meta['Calibration'] = tuple(cals)
        
        data = np.fromfile(self._file_hdl, dtype='<i2', count=1)
        
        # DataType
        meta['DataType'] = data[0]
        
        if not data[0] in self._dictDataType:
            raise RuntimeError('Unknown DataType: "{}"'.format(data[0]))
        if verbose:
            print('DataType:\t{},\t{}'.format(data[0],self._dictDataType[data[0]]))
                
        if self.head['DataTypeID'] == 0x4120:
            # 1D data element
        
            data = np.fromfile(self._file_hdl, dtype='<i4', count=1)
            # ArrayLength
            data = data.tolist()
            meta['ArrayShape'] = data
            if verbose:
                print('ArrayShape:\t{}'.format(data))
                
            dataset = np.fromfile(self._file_hdl, dtype=self._dictDataType[meta['DataType']], count=meta['ArrayShape'][0])
        
        elif self.head['DataTypeID'] == 0x4122:
            # 2D data element
            
            data = np.fromfile(self._file_hdl, dtype='<i4', count=2)
            # ArrayShape
            data = data.tolist()
            meta['ArrayShape'] = data
            if verbose:
                print('ArrayShape:\t{}'.format(data))

            # dataset
            dataset = np.fromfile(self._file_hdl, dtype=self._dictDataType[meta['DataType']], count=meta['ArrayShape'][0]*meta['ArrayShape'][1])
            dataset = dataset.reshape(meta['ArrayShape'][::-1]) #needs to be reversed for little endian data
        
            dataset = np.flipud(dataset)

        return dataset, meta


    def _getTag(self, index, verbose=False):
        '''Retrieve tag from data file.

        Parameters:
            index (int):    Index of tag.
            verbose (bool):    True to get extensive output while reading the file.

        Returns:
            (dict):    Tag as dict.
            
        '''

        # check index, will raise Exceptions if not
        try:
            self._checkIndex(index)
        except:
            raise

        if verbose:
            print('Getting tag {} of {}.'.format(index, self.head['ValidNumberElements']))
        
        tag = {}
            
        try:
            # bad tagoffsets occured pointing to the end of the file
        
            # go to dataset in file
            self._file_hdl.seek(self.head['TagOffsetArray'][index],0)
            
            data = np.fromfile(self._file_hdl, dtype='<i4', count=2)
        
            # TagTypeID
            tag['TagTypeID'] = data[0]
            
            # only proceed if TagTypeID is the same like in the file header (bad TagOffsetArray issue)
            if tag['TagTypeID'] == head['TagTypeID']:
                if verbose:
                    print('TagTypeID:\t"{:#06x}",\t{}'.format(data[0], self._dictTagTypeID[data[0]]))
                    
                # Time    
                tag['Time'] = data[1]
                if verbose:
                    print('Time:\t{}'.format(data[1]))
                
                # check for position
                if tag['TagTypeID'] == 0x4142:
                    data = np.fromfile(self._file_hdl, dtype='<f8', count=2)
                    
                    # PositionX
                    tag['PositionX'] = data[0]
                    if verbose:
                        print('PositionX:\t{}'.format(data[0]))
                    
                    # PositionY
                    tag['PositionY'] = data[1]
                    if verbose:
                        print('PositionY:\t{}'.format(data[1]))
            else:
                # otherwise raise to get to default tag
                raise
                    
        except:
            tag['TagTypeID'] = 0
            tag['Time'] = 0
            tag['PositionX'] = np.nan
            tag['PositionY'] = np.nan
        
     
        return tag
        
    
    def _createDim(self, size, offset, delta, element):
        '''Create dimension labels from SER information.
        
        Parameters:
            size (int):    Number of elements.
            offset (float):    Value at indicated element.
            delta (float):    Difference between elements.
            element (int):    Indicates the element of value offset.
        
        Returns:
            (np.ndarray):    Dimension vector as array.
        
        '''
        
        #dim = np.zeros(size)
        #dim[:] = np.nan
        
        # if element is out off range, map it back into defined
        if element >= size:
            element = size-1
            offset = offset - (element-(size-1))*delta
        
        dim = np.array(range(size)).astype('f8')
        dim = dim*delta
        dim += (offset - dim[element])
        
        # some weird shifting, positionx is +0.5, positiony is -0.5
        # doing this during saving
        # dim += 0.5*delta
        
        return dim
    
    
    def _parseEntry_emi(self, value):
        '''Auxiliary function to parse string entry to int, float or np.string_().
        
        Parameters:
            value (str):    String containing an int, float or string.
        
        Returns:
            (int/float/str):    Value as int, float or string.
            
        '''
        
        # try to parse as int
        try:
            p = int(value)
        except:
            # if not int, then try float
            try:
                p = float(value)
            except:
                # if neither int nor float, stay with string
                p = np.string_(str(value))
        
        return p
        
    
    def read_emi(self, filename):
        '''Read the meta data from an _emi file.
        
        Parameters:
            filename (str):    Name of the _emi file.
        
        Returns:
            (dict):    Dict of metadata.
        '''
        
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string')

        # try opening the file
        try:
            # open file for reading bytes, as binary and text are intermixed
            f_emi = open(filename, 'rb')
        except IOError:
            print('Error reading file: "{}"'.format(filename))
            raise
        except :
            raise
        
        # dict to store _emi stuff
        _emi = {}
        
        # need anything readable from <ObjectInfo> to </ObjectInfo>
        collect = False
        data = b''
        for line in f_emi:
            if b'<ObjectInfo>' in line:
                collect = True
            if collect:
                data += line.strip()
            if b'</ObjectInfo>' in line:
                collect = False

        # close the file
        f_emi.close()
            
        # strip of binary stuff still around
        data = data.decode('ascii', errors='ignore')
        matchObj = re.search('<ObjectInfo>(.+?)</ObjectInfo', data)
        try:
            data = matchObj.group(1)
        except:
            raise RuntimeError('Could not find _emi metadata in specified file.')

        # parse metadata as xml
        root = ET.fromstring('<_emi>'+data+'</_emi>')
        
        # single items
        _emi['Uuid'] = root.findtext('Uuid')
        _emi['AcquireDate'] = root.findtext('AcquireDate')
        _emi['Manufacturer'] = root.findtext('Manufacturer')
        _emi['DetectorPixelHeigth'] = root.findtext('DetectorPixelHeight')
        _emi['DetectorPixelWidth'] = root.findtext('DetectorPixelWidth')

        # Microscope Conditions
        grp = root.find('ExperimentalConditions/MicroscopeConditions')
        
        for elem in grp:
            _emi[elem.tag] = self._parseEntry_emi(elem.text)
            #print('{}:\t{}'.format(elem.tag, elem.text)) 

        # Experimental Description
        grp = root.find('ExperimentalDescription/Root')
        
        for elem in grp:
            _emi['{} [{}]'.format(elem.findtext('Label'), elem.findtext('Unit'))] = self._parseEntry_emi(elem.findtext('Value'))
            #print('{} [{}]: \t{}'.format(elem.findtext('Label'), elem.findtext('Unit'), elem.findtext('Value')))

        # AcquireInfo
        grp = root.find('AcquireInfo')
        
        for elem in grp:
            _emi[elem.tag] = self._parseEntry_emi(elem.text)
            
        # DetectorRange
        grp = root.find('DetectorRange')
        
        for elem in grp:
            _emi['DetectorRange_'+elem.tag] = self._parseEntry_emi(elem.text)

        return _emi
        
        
    def writeEMD(self, filename):
        '''
        Write SER data to an EMD file.
        
        Parameters:
            filename (str):    Name of the EMD file.
        '''
        
        # create the EMD file and set version attributes
        try:
            f = ncempy.io.emd.fileEMD(filename)
        except:
            raise IOError('Cannot write to file "{}"!'.format(filename))
        
        # create EMD group    
        grp = f._file_hdl['data'].create_group(os.path.basename(self._file_hdl.name))
        grp.attrs['emd_group_type'] = 1
            
        
        # use first dataset to layout memory
        data, first_meta = self.getDataset(0)
        first_tag = self._getTag(0)
        
        if self.head['DataTypeID'] == 0x4122:
            # 2D datasets
            
                if first_tag['TagTypeID'] == 0x4142:
                    # 2D mapping
                    dset = grp.create_dataset('data', (self.head['Dimensions'][1]['DimensionSize'], self.head['Dimensions'][0]['DimensionSize'], first_meta['ArrayShape'][1], first_meta['ArrayShape'][0]), dtype=self._dictDataType[first_meta['DataType']] )
                    
                    # collect time
                    time = np.zeros( (self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][1]['DimensionSize']), dtype='i4')

                    # create mapping dims for checking
                    map_xdim = self._createDim( self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][0]['CalibrationOffset'], self.head['Dimensions'][0]['CalibrationDelta'], self.head['Dimensions'][0]['CalibrationElement'] )
                    map_ydim = self._createDim( self.head['Dimensions'][1]['DimensionSize'], self.head['Dimensions'][1]['CalibrationOffset'], self.head['Dimensions'][1]['CalibrationDelta'], self.head['Dimensions'][1]['CalibrationElement'] )
                    # weird direction dependend half pixel shifting
                    map_xdim += 0.5*self.head['Dimensions'][0]['CalibrationDelta']
                    map_ydim -= 0.5*self.head['Dimensions'][1]['CalibrationDelta']

                    for y in range(self.head['Dimensions'][0]['DimensionSize']):
                        for x in range(self.head['Dimensions'][1]['DimensionSize']):

                            index = int(x+y*self.head['Dimensions'][0]['DimensionSize'])
                            print('converting dataset {} of {}, items ({}, {})'.format(index+1, self.head['ValidNumberElements'], x, y))
                   
                            # retrieve dataset and put into buffer     
                            data, meta = self.getDataset( index )
                            dset[y, x, :,:] = data[:,:]
                            
                            # get tag data per image
                            tag = self._getTag(index)
                            time[y,x] = tag['Time']

                            assert( np.abs(tag['PositionX'] - map_xdim[x]) < np.abs(tag['PositionX']*1e-8) )
                            assert( np.abs(tag['PositionY'] - map_ydim[y]) < np.abs(tag['PositionY']*1e-8) )
                    
                            del data, meta, tag
                    
                    # create dimension datasets
                    dims = []
                    dims_time = []
                    
                    # Position Y
                    assert self.head['Dimensions'][1]['Description'] == 'Position'
                    dims.append( (map_ydim, self.head['Dimensions'][1]['Description'], '[{}]'.format(self.head['Dimensions'][1]['Units'])) )
                    dims_time.append( (map_ydim, self.head['Dimensions'][1]['Description'], '[{}]'.format(self.head['Dimensions'][1]['Units'])) ) 
                    
                    # Position X
                    assert self.head['Dimensions'][0]['Description'] == 'Position'
                    dims.append( (map_xdim, self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )
                    dims_time.append( (map_xdim, self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )
                    
                    dim = self._createDim(first_meta['ArrayShape'][1], first_meta['Calibration'][1]['CalibrationOffset'], first_meta['Calibration'][1]['CalibrationDelta'], first_meta['Calibration'][1]['CalibrationElement'])
                    dims.append( (dim, 'y', '[m]') )
                    
                    dim = self._createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
                    dims.append( (dim, 'x', '[m]') )
                    
                    # write dimensions
                    for i in range(len(dims)):
                        f.write_dim('dim{:d}'.format(i+1), dims[i], grp)
 
                    # write out time as additional dataset
                    grp = f.put_emdgroup('timestamp', time, dims_time, parent=grp)
                    
                    
                else:
                
                    # 1 entry series to single image
                    if self.head['ValidNumberElements'] == 1:
                    
                        # get image
                        data, meta = self.getDataset(0)
                        tag = self._getTag(0)
                        
                        # create dimensions
                        dims = []
                        
                        dim = self._createDim(first_meta['ArrayShape'][1], first_meta['Calibration'][1]['CalibrationOffset'], first_meta['Calibration'][1]['CalibrationDelta'], first_meta['Calibration'][1]['CalibrationElement'])
                        dims.append( (dim, 'y', '[m]') )
                                                                          
                        dim = self._createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
                        dims.append( (dim, 'x', '[m]') )
                        
                        dset = grp.create_dataset( 'data', (first_meta['ArrayShape'][1], first_meta['ArrayShape'][0]), dtype=self._dictDataType[first_meta['DataType']])
                        
                        dset[:,:] = data[:,:]
                        
                        for i in range(len(dims)):
                            f.write_dim('dim{:d}'.format(i+1), dims[i], grp)
                            
                        dset.attrs['timestamp'] = tag['Time']
                        
                    else:
                
                        # simple series
                        dset = grp.create_dataset( 'data', (self.head['ValidNumberElements'], first_meta['ArrayShape'][1], first_meta['ArrayShape'][0]), dtype=self._dictDataType[first_meta['DataType']] )
            
                        # collect time
                        time = np.zeros(self.head['ValidNumberElements'], dtype='i4')
            
                        for i in range(self.head['ValidNumberElements']):
                            print('converting dataset {} of {}'.format(i+1, self.head['ValidNumberElements']))
                
                            # retrieve dataset and put into buffer
                            data, meta = self.getDataset(i)
                            dset[i,:,:] = data[:,:]
                
                            # get tag data per image
                            tag = self._getTag(i)
                            time[i] = tag['Time']
                            
                        # create dimension datasets
                        dims = []
                            
                        # first SER dimension is number
                        assert self.head['Dimensions'][0]['Description'] == 'Number'
                            
                        dim = self._createDim(self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][0]['CalibrationOffset'], self.head['Dimensions'][0]['CalibrationDelta'], self.head['Dimensions'][0]['CalibrationElement'])
                        dims.append( (dim[0:self.head['ValidNumberElements']], self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )

                        dim = self._createDim(first_meta['ArrayShape'][1], first_meta['Calibration'][1]['CalibrationOffset'], first_meta['Calibration'][1]['CalibrationDelta'], first_meta['Calibration'][1]['CalibrationElement'])
                        dims.append( (dim, 'y', '[m]') )
                                                                          
                        dim = self._createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
                        dims.append( (dim, 'x', '[m]') )

                        # write dimensions
                        for i in range(len(dims)):
                            f.write_dim('dim{:d}'.format(i+1), dims[i], grp)
                        
                        # write out time as additional dim vector
                        f.write_dim('dim1_time', (time, 'timestamp', '[s]'), grp)
        
        
        elif self.head['DataTypeID'] == 0x4120:
            # 1D datasets
            
            if first_tag['TagTypeID'] == 0x4142:
                    # 2D mapping
                    dset = grp.create_dataset( 'data', (self.head['Dimensions'][1]['DimensionSize'], self.head['Dimensions'][0]['DimensionSize'], first_meta['ArrayShape'][0]), dtype=self._dictDataType[first_meta['DataType']] )
                    
                    time = np.zeros( (self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][1]['DimensionSize']), dtype='i4')

                    # create mapping dims for checking
                    map_xdim = self._createDim( self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][0]['CalibrationOffset'], self.head['Dimensions'][0]['CalibrationDelta'], self.head['Dimensions'][0]['CalibrationElement'] )
                    map_ydim = self._createDim( self.head['Dimensions'][1]['DimensionSize'], self.head['Dimensions'][1]['CalibrationOffset'], self.head['Dimensions'][1]['CalibrationDelta'], self.head['Dimensions'][1]['CalibrationElement'] )
                    # weird direction dependend half pixel shifting
                    map_xdim += 0.5*self.head['Dimensions'][0]['CalibrationDelta']
                    map_ydim -= 0.5*self.head['Dimensions'][1]['CalibrationDelta']

                    for y in range(self.head['Dimensions'][0]['DimensionSize']):
                        for x in range(self.head['Dimensions'][1]['DimensionSize']):

                            index = int(x+y*self.head['Dimensions'][0]['DimensionSize'])
                            print('converting dataset {} of {}, items ({}, {})'.format(index+1, self.head['ValidNumberElements'], x, y))
                   
                            # retrieve dataset and put into buffer     
                            data, meta = self.getDataset( index )
                            dset[y, x, :] = np.copy(data[:])
                            
                            # get tag data per image
                            tag = self._getTag(index)
                            time[y,x] = tag['Time']
                                
                            assert( np.abs(tag['PositionX'] - map_xdim[x]) < np.abs(tag['PositionX']*1e-8) )
                            assert( np.abs(tag['PositionY'] - map_ydim[y]) < np.abs(tag['PositionY']*1e-8) )
                    
                            del data, meta, tag
                    
                    # create dimension datasets
                    dims = []
                    dims_time = []
                    
                    # Position Y
                    assert self.head['Dimensions'][1]['Description'] == 'Position'
                    dims.append( (map_ydim, self.head['Dimensions'][1]['Description'], '[{}]'.format(self.head['Dimensions'][1]['Units'])) )
                    dims_time.append( (map_ydim, self.head['Dimensions'][1]['Description'], '[{}]'.format(self.head['Dimensions'][1]['Units'])) )
                    
                    # Position X
                    assert self.head['Dimensions'][0]['Description'] == 'Position'
                    dims.append( (map_xdim, self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )
                    dims_time.append( (map_xdim, self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )
                    
                    dim = self._createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
                    dims.append( (dim, 'E', '[m_eV]') )
                        
                    # write dimensions
                    for i in range(len(dims)):
                        f.write_dim('dim{:d}'.format(i+1), dims[i], grp)
 
                    # write out time as additional dataset
                    grp = f.put_emdgroup('timestamp', time, dims_time, parent=grp)
                    
            else:
                    # simple series
                    dset = grp.create_dataset( 'data', (self.head['ValidNumberElements'], first_meta['ArrayShape'][0]), dtype=self._dictDataType[first_meta['DataType']])
                    
                    # collect time
                    time = np.zeros(self.head['ValidNumberElements'], dtype='i4')
                    
                    for i in range(self.head['ValidNumberElements']):
                        print('converting dataset {} of {}'.format(i+1, self.head['ValidNumberElements']))
            
                        # retrieve dataset and put into buffer
                        data, meta = self.getDataset(i)
                        dset[i,:] = data[:]
            
                        # get tag data per image
                        tag = self._getTag(i)
                        time[i] = tag['Time']

                    # create dimension datasets
                    dims = []

                    # first SER dimension is number
                    assert self.head['Dimensions'][0]['Description'] == 'Number'
                    dim = self._createDim(self.head['Dimensions'][0]['DimensionSize'], self.head['Dimensions'][0]['CalibrationOffset'], self.head['Dimensions'][0]['CalibrationDelta'], self.head['Dimensions'][0]['CalibrationElement'])
                    dims.append( (dim[0:self.head['ValidNumberElements']], self.head['Dimensions'][0]['Description'], '[{}]'.format(self.head['Dimensions'][0]['Units'])) )
                        
                    dim = self._createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
                    dims.append( (dim, 'E', '[m_eV]') )
                        
                    # write dimensions
                    for i in range(len(dims)):
                        f.write_dim('dim{:d}'.format(i+1), dims[i], grp)
 
                    # write out time as additional dim vector
                    f.write_dim('dim1_time', (time, 'timestamp', '[s]'), grp)    
                    
        else:
            raise RuntimeError('Unknown DataTypeID')    
            

        # attempt to free memory asap
        #del dset_buf, dims
        
        
        # put meta information from _emi to Microscope group, if available
        if self._emi:
            for key in self._emi:
                if not self._emi[key] is None:
                    f.microscope.attrs[key] = self._emi[key]
                
        # write comment into Comment group
        f.put_comment('Converted SER file "{}" to EMD using the openNCEM tools.'.format(self._file_hdl.name))
        
def serReader(fName):
    '''Simple function to parse the file and read all datasets. This is a one function implementation to load all data in a ser file.
    '''
    f1 = fileSER(fName) #open the file and init the class
    if f1.head['ValidNumberElements'] > 0:
        im1 = f1.getDataset(0) #im1 will be used as the output dictionary to keep most of the meta data
        im1[1]['filename'] = fName #save the file name in the output dictionary
        temp = np.empty([f1.head['ValidNumberElements'], im1[0].shape[0],im1[0].shape[1]],dtype=im1[0].dtype)
        temp[0,:,:] = im1[0]
        for ii in range(1,f1.head['ValidNumberElements']):
            temp[ii,:,:] = f1.getDataset(ii)[0] #get the next dataset
        im1[1]['data'] = temp
    else:
        im1 = {}
        print('No data set found') 
    del f1 #delete the class and close the file
    return im1[1] #return the dataset and metadata as a dictionary