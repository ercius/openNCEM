"""
A module to read MRC files in python and numpy.
Written according to MRC specification at http://bio3d.colorado.edu/imod/betaDoc/mrc_format.txt
Also works with FEI MRC files which include a special header block with experimental information.

Note
----
General users:
    Use the simplified mrc.mrcReader() function to load the data and meta
    data as a python dictionary.

Advanced users and developers:
    Access the file internals through the mrc.fileMRC() class.

"""

from pathlib import Path

import numpy as np


class fileMRC:
    """ Read in the data in MRC format and other useful information like metadata. Follows the specification
    published at http://bio3d.colorado.edu/imod/betaDoc/mrc_format.txt

    Attributes
    ----------
    file_name : str
        The name of the file
    file_path : pathlib.Path
        A pathlib.Path object for the open file
    fid : file
        The file handle to the opened MRC file.
    mrcType : int
        The internal MRC data type.
    dataType : np.dtype
        The numpy dtype corresponding to the mrcType.
    dataSize : np.ndarray
        The number of pixels along each dimension. Corresponds to the shape attribute of a np.ndarray
    gridSize : np.ndarray
        The size of the grid. Usually the same as dataSize
    volumeSize : np.ndarray
        The size of the volume along each direction in Angstroms.
    voxelSize : np.ndarray
        The size of the voxel along each direction in Angstroms.
    cellAngles : np.ndarray
        The angles of the cell. Ignored in most cases including ncempy.
    axisOrientations : np.ndarray
        Mapping the orientations of the data to real space directions X, Y, Z. Ignored by ncempy
    minMaxMean  : np.ndarray
        The minimum, maximum and mean value of the data to avoid computing this every time.
    extra : np.ndarray
        Extra mbinary metadata if it exists.
    FEIinfo : dict
        A dictionary of metadata used by FEI (Thermo Fischer) microsocpes for important metadata. This metadata
        overwrites the voxelsize attribute if it exists.
    dataOffset : int
        The integer offset in bytes to the start of the raw data.
    dataOut: dict
        Will hold the data and metadata to output to the user after getDataset() call.
    v : bool
        More output for debugging. False by default

    Examples
    --------
    Read in all data and metadata into memory.
    >> import ncempy.io as nio
    >> mrc0 = nio.mrc.mrcReader('file.mrc')

    Low level operations to get 1 slice of the 3D data
    >> import ncempy.io as nio
    >> with nio.mrc.fileMRC('file.mrc') as f1:
    >>     single_slice = f1.getSlice(0)
    """

    def __init__(self, filename, verbose=False):
        """
        Parameters
        -----------
            filename : str or pathlib.Path or file object
                String or pathlib.Path of file object pointing to the filesystem location of the file.
            verbose : bool
                If True, debug information is printed.

        """
        # check filename type
        if hasattr(filename, 'read'):
            self.fid = filename
            try:
                self.file_path = Path(self.fid.name)
                self.file_name = self.file_path.name
            except AttributeError:
                self.file_name = None
                self.file_path = None
        else:
            if isinstance(filename, str):
                filename = Path(filename)
            elif isinstance(filename, Path):
                pass
            else:
                raise TypeError('Filename is supposed to be a string or pathlib.Path or readable file object')

            self.file_path = filename
            self.file_name = filename.name

            # Open the file and quit if the file does not exist
            try:
                self.fid = open(self.file_path, 'rb')
            except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))

        # necessary declarations, if something fails
        self.mrcType = None
        self.dataType = None
        self.dataSize = None
        self.gridSize = None
        self.volumeSize = None
        self.voxelSize = None
        self.cellAngles = None
        self.axisOrientations = None
        self.minMaxMean = None
        self.extra = None
        self.FEIinfo = None
        self.dataOffset = None

        self.dataOut = {}  # will hold the data and metadata to output to the user after getDataset() call

        # Add a top level variable to indicate verbose output for debugging
        self.v = verbose

        self.parseHeader()

    def __del__(self):
        """Close the file.

        """
        if not self.fid.closed:
            if self.v:
                print('Closing input file: {}'.format(str(self.file_path)))
            self.fid.close()
        return None

    def __enter__(self):
        """Implement python's with statement

        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Implement python's with statement.
        Close the file using __del__

        """
        self.__del__()
        return None

    def parseHeader(self):
        """Read the header information which includes data type, data size, data
        shape, and metadata.

        Note
        -----
            This header uses Fortran-style ordering. Numpy uses C-style ordering.
            The header is read in and then some attributes are reversed [::-1] at
            the end for output to the user to match C-ordering in numpy.

        TODO
        ----
            Implement special dtype to read the entire header at once.
            Read everything at once using special dtype. ~5x faster than multiple np.fromfile() reads:

            Untested but works in theory
            headerDtype = np.dtype([('head1','10int32'),('head2','6float32'),('axisOrientations','3int32'),
            ('minMaxMean','3int32'),('extra','32int32')])
            head = np.fromfile(self.fid,dtype=headerDtype,count=1)
        """

        # Always start at the beginning of the file.
        self.fid.seek(0)

        header_dtype = np.dtype(
            [('head1', '10int32'), ('head2', '6float32'), ('axisOrientations', '3int32'), ('minMaxMean', '3int32'),
             ('extra', '32int32')])
        header = np.fromfile(self.fid, dtype=header_dtype, count=1)

        # Read in the initial header values
        head1 = header['head1'][0]
        if self.v:
            print('header1 = {}'.format(head1))
        # Set the number of pixels for each dimension
        self.dataSize = head1[0:3]
        if self.v:
            print('dataSize (fortran ordering) = {}'.format(self.dataSize))

        # Set the data type and convert to numpy type
        self.mrcType = head1[3]
        self.dataType = self._getMRCType(self.mrcType)
        if self.v:
            print('dataType = {}'.format(self.dataType))

        # Get the grid size
        self.gridSize = head1[7:10]
        if self.v:
            print('mrc defined gridSize = {}'.format(self.gridSize))

        # Get the physical volume size (always in Angstroms) (starting at byte #11 in the file).
        head2 = header['head2'][0]
        self.volumeSize = head2[0:3]
        if self.v:
            print('mrc defined volumeSize = {}'.format(self.volumeSize))

        # calculate the voxel size based on volume and grid sizes
        # account for zero values which are often stored in volume and grid size values
        self.voxelSize = np.ones(3)
        for ii in range(0, 3):
            if self.volumeSize[ii] == 0:
                v = 1
            else:
                v = self.volumeSize[ii]
            if self.gridSize[ii] == 0:
                gs = 1
            else:
                gs = self.gridSize[ii]
            self.voxelSize[ii] = v / gs

        """if self.volumeSize.any() == 0 or self.gridSize.any() == 0:
            if self.v:
                print('Detected 0 volume or grid size. Setting voxel size to 1 (Ang).')
            self.voxelSize = np.ones(3)  # use 1 as a voxel size if its not set in the file
        else:
            self.voxelSize = (self.volumeSize / np.float32(self.gridSize))
            if self.v:
                print('voxelSize (Ang) = {}'.format(self.voxelSize))
        """

        # Pixel (cell) angles
        self.cellAngles = head2[3:6]
        if self.v:
            print('cellAngles = {}'.format(self.cellAngles))

        # Axis orientations. Tells which axes are X,Y,Z
        self.axisOrientations = header['axisOrientations'][0]
        if self.v:
            print('axisOrientations = {}'.format(self.axisOrientations))

        # Min, max,mean
        self.minMaxMean = header['minMaxMean'][0]

        # Extra information (for FEI MRC file, extra(1) is the size of the FEI information encoded with the file in
        # terms of 4 byte floats)
        self.extra = header['extra'][0]

        # Numpy uses C-style ordering. The header is written in Fortran-Style ordering.
        # Flip the order of everything useful
        if self.v:
            print(
                'Note: The MRC header is written in Fortran-Style ordering, but Numpy uses C-style ordering. '
                'This program will now reverse the order (using [::-1]) of useful metadata: (dataSize, gridSize, '
                'volumeSize, voxelSize, cellAngles, axisOrientations)')

        self.dataSize = self.dataSize[::-1]
        self.gridSize = self.gridSize[::-1]
        self.volumeSize = self.volumeSize[::-1]
        self.voxelSize = self.voxelSize[::-1]
        self.cellAngles = self.cellAngles[::-1]
        self.axisOrientations = self.axisOrientations[::-1]
        # self.Shape = self.Shape[::-1]

        # Move to the end of the normal header
        self.fid.seek(1024)

        # Read in the extended header if it exists (for FEI MRC files)
        if self.extra[1] != 0:
            pos1 = self.fid.tell()
            if self.v:
                print('Extra header found. Most likely and FEI-style MRC file.')
                print('Position before reading extra header: ' + str(pos1))
                print('Extra header size = ' + str(self.extra[1]))

            ''' Read the extra FEI header described as follows:
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
            FEIinfoValues = np.fromfile(self.fid, dtype=np.float32, count=15)
            self.FEIinfo = {'a_tilt': FEIinfoValues[0], 'b_tilt': FEIinfoValues[1], 'x_stage': FEIinfoValues[2],
                            'y_stage': FEIinfoValues[3], 'z_stage': FEIinfoValues[4], 'x_shift': FEIinfoValues[5],
                            'y_shift': FEIinfoValues[6], 'defocus': FEIinfoValues[7], 'exposure_time': FEIinfoValues[8],
                            'mean': FEIinfoValues[9], 'tilt_axis': FEIinfoValues[10], 'pixel_size': FEIinfoValues[11],
                            'magnification': FEIinfoValues[12], 'voltage': FEIinfoValues[13],
                            'unknown': FEIinfoValues[14]}

            self.voxelSize[0] = 1.  # set this to 1 but it should be the tilt angles. These can be non-uniform though.
            if self.FEIinfo['pixel_size'] != 0:
                self.voxelSize[1] = self.FEIinfo['pixel_size'] * 1e10  # convert meter to Angstroms, standard for MRCs
                self.voxelSize[2] = self.FEIinfo['pixel_size'] * 1e10

            if self.v:
                print('Extended header data')
                for aa, bb in self.FEIinfo.items():
                    print('{} = {}'.format(aa, bb))
        else:
            self.FEIinfo = {}

        self.dataOffset = 1024 + self.extra[1]  # offset of the data from the start of the file

        # Add relevant information (metadata) to the output dictionary
        self.dataOut = {'pixelSize': self.voxelSize, 'voxelSize': self.voxelSize,
                        'cellAngles': self.cellAngles, 'axisOrientations': self.axisOrientations,
                        'filename': self.file_name}
        if self.extra[1] != 0:
            self.dataOut['FEIinfo'] = self.FEIinfo

        return 1

    def getDataset(self):
        """Read in the full data block and reshape to an ndarray
        with C-style ordering.

        """
        self.fid.seek(self.dataOffset, 0)  # move to the start of the data from the start of the file
        try:
            num0 = int(np.prod(self.dataSize, dtype=np.uint64))
            data1 = np.fromfile(self.fid, dtype=self.dataType, count=num0)
            self.dataOut['data'] = data1.reshape(self.dataSize)
        except MemoryError:
            print("Not enough memory to read in the full data set. Use getMemmap")
        return self.dataOut

    def getSlice(self, num):
        """Read in a slice of an MRC file. Useful for parsing through a large file without reading
        the entire data set into memory.

        Parameters
        ----------
            num : int
                Get the requested image.

        Returns
        -------
            out : ndarray
                A 2D slice or a 3D set of slices along the first index

        Raises
        ------
            IndexError
                If num > the number of slices.

        """
        # Check num is within the data array size bounds
        if num > (self.dataSize[0] - 1):
            raise IndexError('Index {} is out of bounds for array with size {}'.format(num, self.dataSize[0]))

        self.fid.seek(self.dataOffset, 0)  # move to the start of the data by skipping the header
        imSize = self.dataSize[1] * self.dataSize[2]  # size of each image in pixels
        byteSize = int(np.dtype(self.dataType).itemsize * imSize)  # make sure this is an int32
        self.fid.seek(num * byteSize, 1)  # skip to the slice requested from the start of the data

        data1 = np.fromfile(self.fid, dtype=self.dataType, count=imSize)  # read in the requested image
        data1 = data1.reshape((self.dataSize[1], self.dataSize[2]))  # reshape the image

        return data1

    def getMemmap(self):
        """Return a numpy memmap object (read-only) for the dataset. This is very useful
        for very large datasets to avoid loading the entire data set into memory. No meta data is
        returned.

        Returns
        --------
        : numpy.core.memmap
            A read-only numpy memmap object with access to the data on disk.
        """
        mm = np.memmap(self.fid, dtype=self.dataType, mode='r', offset=self.dataOffset,
                       shape=tuple(self.dataSize))

        return mm

    def _applyAxisOrientations(self, arrayIn):
        """ This is untested and unused.

        """
        return [arrayIn[x - 1] for x in self.axisOrientations]

    def _getMRCType(self, dataType):
        """Return the correct data type according to the official MRC type list:
        0 image : signed 8-bit bytes range -128 to 127
        1 image : 16-bit
        2 image : 32-bit reals
        3 transform : complex 16-bit integers
        4 transform : complex 32-bit reals
        6 image : unsigned 16-bit range 0 to 65535

        Parameters
        ----------
            dataType: int
                The data type value encoded in an MRC header

        Returns
        --------
            out: numpy dtype
                The corresponding numpy data type.

        """
        if dataType == 0:
            Type = np.int8
        elif dataType == 1:
            Type = np.int16
        elif dataType == 2:
            Type = np.float32
        elif dataType == 6:
            Type = np.uint16
        else:
            print("Unsupported data type" + str(dataType))  # complex data types are currently unsupported
            Type = None
        return Type


def mrcReader(file_name):
    """A simple function to read open a MRC, parse the header, and read the full
    data set.

    Parameters
    ----------
        file_name : str or pathlib.Path
            The path to the file to load.

    Returns
    -------
        out: dict
            A dictionary containing the data and interesting metadata. The data is attached to the 'data' key.

    Example
    -------
        Read in all data from disk into memory. This assumes the dataset is 3 dimensional:
        >> from ncempy.io.mrc import mrcReader
        >> import matplotlib.pyplot as plt
        >> mrc1 = mrcReader('filename.mrc')
        >> plt.imshow(mrc1['data'][0, :, :])  # show the first image in the data set
    """
    if isinstance(file_name, str):
        file_name = Path(file_name)

    with fileMRC(file_name) as f1:  # open the file and init the class
        im1 = f1.getDataset()  # read in the dataset
    
    # Add extra meta data not already in im1
    im1['filename'] = file_name.name
    im1['pixelUnit'] = 'A'

    return im1  # return the data and metadata as a dictionary


def mrc2raw(file_name):
    """ Convert an MRC type file to raw binary. The entire data is read into memory and then written to a new
    raw file in the same location with the data type and shape (C-ordering) written into the filename. No other
    meta data is retained.

    Parameters
    ----------
    file_name : str or pathlib.Path
        The file name to convert.

    """
    if isinstance(file_name, str):
        file_name = Path(file_name)

    tomo = mrcReader(file_name)

    # Add metadata to file name
    out_name = file_name.stem + '_{}_{}.raw'.format(tomo['data'].dtype, tomo['data'].shape)

    outPath = file_name.with_name(out_name.replace(' ', ''))
    print(outPath)

    with open(outPath, 'wb') as f0:
        f0.write(tomo['data'])  # write out as C ordered data


def mrc2emd(file_name):
    """Write an MRC file as an HDF5 file in EMD format with same file name and .emd ending.
    Header information is retained as attributes.

    Parameters
    ----------
        file_name: str
            The name of the file to convert from MRC to EMD format.

    Returns
    -------
        out: int
            1 if successful.

    TODO
    -----
        Update this to use ncempy.emd class
    """
    import h5py

    # Read in the MRC data and reshape to C-style ordering
    tomo = mrcReader(file_name)

    # create the HDF5 file
    with h5py.File(file_name.rsplit('.mrc', 1)[0] + '.emd', 'w') as f1:  # w- will error if the file exists

        # Create the axis vectors in nanometers. Standard MRC pixel size is in Angstroms
        xFull = np.linspace(0, tomo['voxelSize'][0] * tomo['data'].shape[0] - 1, tomo['data'].shape[0])
        yFull = np.linspace(0, tomo['voxelSize'][1] * tomo['data'].shape[1] - 1, tomo['data'].shape[1])
        zFull = np.linspace(0, tomo['voxelSize'][2] * tomo['data'].shape[2] - 1, tomo['data'].shape[2])

        # Root data group
        dataTop = f1.create_group('data')

        # Create tilt series group
        tiltseriesGroup = dataTop.create_group('data')
        tiltseriesGroup.attrs['emd_group_type'] = np.int8(1)

        # Save the data to the EMD file and reshape it to a C-style array
        try:
            _ = tiltseriesGroup.create_dataset('data', data=tomo['data'], compression='gzip', shuffle=True)
        except MemoryError:
            raise MemoryError("Not enough memory to write out data to EMD file")

        dim1 = tiltseriesGroup.create_dataset('dim1', data=xFull)
        dim1.attrs['name'] = np.string_('x')
        dim1.attrs['units'] = np.string_('')
        dim2 = tiltseriesGroup.create_dataset('dim2', data=yFull)
        dim2.attrs['name'] = np.string_('y')
        dim2.attrs['units'] = np.string_('')
        dim3 = tiltseriesGroup.create_dataset('dim3', data=zFull)
        dim3.attrs['name'] = np.string_('z')
        dim3.attrs['units'] = np.string_('')

        # Create the other groups
        scopeGroup = f1.create_group('Microscope')
        scopeGroup.attrs['voxel sizes'] = tomo['voxelSize']
        userGroup = f1.create_group('User')
        commentGroup = f1.create_group('Comments')


def mrcWriter(filename, data, pixelSize, forceWrite=False):
    """Write out a MRC type file according to the specification at http://bio3d.colorado.edu/imod/doc/mrc_format.txt

    Parameters
    ----------
        filename : str or pathlib.Path
            The name or Path of the file to write out to.
        data : ndarray
            The array data to write to disk.
        pixelSize : tuple
            The size of the pixel along each direction (in Angstroms) as a 3 element vector (sizeZ,sizeY,sizeX).
        forceWrite : bool
            This will write the data as a C-contiguous array. It is not suggested to use this option and it will
            be removed in future versions.

    """

    with open(filename, 'wb') as fid:

        if len(data.shape) > 3:
            print("Too many dimensions")
            return

        if not data.flags['C_CONTIGUOUS']:
            print("Error: Array must be C-style ordering: [numImages,Y,X]. "
                  "Use numpy.tranpspose and np.ascontiguousarray to change data ordering in memory")
            print('Exiting')
            return

        # initialize the header with 256 zeros with size 4 bytes
        header = np.zeros(256, dtype=np.int32)
        fid.write(header)
        fid.seek(0, 0)  # return to the beginning of the file

        # Initialize the int32 part of the header
        header1 = np.zeros(10, dtype=np.int32)

        # Write the number of columns, rows and sections (images)
        # header1[0:3] = np.int32(dims) #data size in pixels
        header1[0] = np.int32(data.shape[2])  # num columns, the last index in C-style ordering
        header1[1] = np.int32(data.shape[1])  # num rows
        header1[2] = np.int32(data.shape[0])  # num sections (images)

        if data.dtype == np.float32:
            header1[3] = np.int32(2)
        elif data.dtype == np.uint16:
            header1[3] = np.int32(6)
        elif data.dtype == np.int16:
            header1[3] = np.int32(1)
        elif data.dtype == np.int8:
            header1[3] = np.int32(0)
        else:
            print("Data type {} is unsupported. Only int8, int16, uint16, and float32 are supported".format(data.dtype))
            return

        # Starting point of sub image (not used in IMOD)
        header1[4:7] = np.zeros(3, dtype=np.int32)

        # Grid size in X,Y,Z
        # header1[7:10] = np.int32(dims); #data size in pixels
        header1[7] = np.int32(data.shape[2])  # mx
        header1[8] = np.int32(data.shape[1])  # my
        header1[9] = np.int32(data.shape[0])  # mz

        # Write out the first part of the header information
        fid.write(header1)

        # Cell dimensions (in Angstroms)
        # pixel spacing = xlen/mx, ylen/my, zlen/mz
        fid.write(np.float32(pixelSize[2] * data.shape[2]))  # xlen
        fid.write(np.float32(pixelSize[1] * data.shape[1]))  # ylen
        fid.write(np.float32(pixelSize[0] * data.shape[0]))  # zlen

        # Cell angles (in degrees)
        fid.write(np.float32([90.0, 90.0, 90.0]))

        # Description of array directions with respect to: Columns, Rows, Images
        fid.write(np.int32([1, 2, 3]))

        # Minimum and maximum density
        fid.write(np.float32(np.min(data)))
        fid.write(np.float32(np.max(data)))
        fid.write(np.float32(np.mean(data)))

        # Needed to indicate that the data is little endian for NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above
        fid.seek(212, 0)
        fid.write(np.int8([68, 65, 0, 0]))  # use [17,17,0,0] for big endian

        # Write out the data
        fid.seek(1024)
        if forceWrite:
            fid.write(np.ascontiguousarray(data))  # Change to C ordering array for writing to disk
        else:
            fid.write(data)  # must be C-contiguous

        # Close the file


def writeHeader(filename, shape, dtype, pixelSize):
    """Write out a MRC type file header according to the specification at
    http://bio3d.colorado.edu/imod/doc/mrc_format.txt. This is useful for initializing an MRC file and then writing to
    it manually or see appendData() function below.

    Parameters
    ----------
    filename : str
        The name of the EMD file
    shape : tuple
        The shape of the data to write
    dtype : numpy.dtype
        The dtype to write out the data as. Only some numpy dtypes are supported byt his format. It is suggested
        to use np.float32 in most cases for maximum compatibility.
    pixelSize : tuple
        The size of the pixel along each direction (in Angstroms) as a 3 element vector (sizeX,sizeY,sizeZ).
        sizeZ could be the angular step for a tilt series

    """

    with open(filename, 'wb') as fid:
        if len(shape) > 3:
            print("Too many dimensions")
            return 0

        # Initialize the header with 256 zeros with size 4 bytes
        header = np.zeros(256, dtype=np.int32)
        fid.write(header)
        fid.seek(0, 0)  # return to the beginning of the file

        # Initialize the int32 part of the header with zeros
        header1 = np.zeros(10, dtype=np.int32)

        # Write the number of columns, rows and sections (images)
        header1[0] = np.int32(shape[2])  # num columns, the last index in C-style ordering
        header1[1] = np.int32(shape[1])  # num rows
        header1[2] = np.int32(shape[0])  # num sections (images)

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
            return 0

        # Starting point of sub image (not used in IMOD)
        header1[4:7] = np.zeros(3, dtype=np.int32)

        # Grid size in X,Y,Z
        header1[7] = np.int32(shape[2])  # mx
        header1[8] = np.int32(shape[1])  # my
        header1[9] = np.int32(shape[0])  # mz

        # Write out the first part of the header information
        fid.write(header1)

        # Cell dimensions (in Angstroms)
        # pixel spacing = xlen/mx, ylen/my, zlen/mz
        fid.write(np.float32(pixelSize[2] * shape[2]))  # xlen
        fid.write(np.float32(pixelSize[1] * shape[1]))  # ylen
        fid.write(np.float32(pixelSize[0] * shape[0]))  # zlen

        # Cell angles (in degrees)
        fid.write(np.float32([90.0, 90.0, 90.0]))

        # Description of array directions with respect to: Columns, Rows, Images
        fid.write(np.int32([1, 2, 3]))

        # Minimum and maximum density
        # fid.write(np.float32(np.min(data)))
        # fid.write(np.float32(np.max(data)))
        # fid.write(np.float32(np.mean(data)))
        fid.write(np.zeros(3, dtype=np.float32))

        # Needed to indicate that the data is little endian for NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above
        fid.seek(212, 0)
        fid.write(np.int8([68, 65, 0, 0]))  # use [17,17,0,0] for big endian


def appendData(filename, data):
    """Append a binary set of data to the end of a MRC file. This should only be used in conjunction with
    writeHeader() above.

    Parameters
    ----------
    filename : str
        Name of the MRC file with pre-initiated header and some data already written.
    data : ndarray
        Data to append to the file.

    """
    with open(filename, 'ab') as fid:
        # Write out the data
        fid.seek(0, 2)  # seek to the end of the file
        fid.write(data)  # Change to C ordering array for writing to disk


def emd2mrc(filename, dsetPath):
    """Convert EMD data set into MRC data set. The final data type is float32 for convenience.

    Parameters
    ----------
    filename : str
        The name of the EMD file
    dsetPath : str
        The HDF5 path to the top group holding the data. ex. '/data/raw/'
    """
    import h5py
    with h5py.File(filename, 'r') as f1:
        # Get the pixel sizes and convert to Ang
        dimsPath = dsetPath + '/dim'
        print('Warning: Assuming EMD dim vectors are in nanometer')
        print('dim vector names/units are:')
        for ii in range(1, 4):
            print('name, units = {}, {}'.format(f1[dimsPath + str(ii)].attrs['name'],
                                                f1[dimsPath + str(ii)].attrs['units']))
        pixelSizeX = (f1[dsetPath + '/dim2'][1] - f1[dsetPath + '/dim2'][0]) * 10  # change nanometers to Ang
        pixelSizeY = (f1[dsetPath + '/dim3'][1] - f1[dsetPath + '/dim3'][0]) * 10  # change nanometers to Ang

        filenameOut = filename.split('.emd')[
                          0] + '.mrc'  # use the first part of the file as the prefix removing the .emd on the end

        print('Warning: Converting to float32 before writing to disk')
        mrcWriter(filenameOut, np.float32(f1[dsetPath + '/data']),
                  (1, pixelSizeY, pixelSizeX))

        print('Finished writing to: {}'.format(filenameOut))
