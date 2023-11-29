"""
This module provides an interface to the Berkeley EMD file format.

See https://emdatasets.com/ for more details.

* It fully supports v0.2 for reading and writing

Note
----
General users:
    Use the simplified emd.emdReader() function to load the data and meta
    data as a python dictionary.

Advanced users and developers:
    Access the file internals through the emd.fileDEMD() class.

"""

from pathlib import Path
import datetime

import numpy as np
import h5py


class NoEmdDataSets(Exception):
    """Special exception indicating no EMD datasets are in an EMD file."""
    pass

class UnsupportedEmdVersion(Exception):
    """Special exception indicating an unsupported EMD version."""
    pass
        


class fileEMD:
    """Class to represent Berkeley EMD files.

    Implemented for spec 0.2 using the recommended layout for metadata.

    Attributes
    ----------
    file_name : str
        The name of the file
    file_path : pathlib.Path
        A pathlib.Path object for the open file
    file_hdl : h5py.File
        The h5py file handle which provides direct access to the underlying h5py file structure.
    version : tuple
        A 2-tuple of the major and minor EMD version as ints.
    data : h5py.Group
        A link to the data group in the EMD file.
    microscope : h5py.Group
        A link to the microscope EMD group with metadata.
    sample : h5py.Group
        A link to the sample EMD group with metadata.
    user : h5py.Group
        A link to the user EMD group with metadata.
    comments : h5py.Group
        A link to the comments in this EMD file. A time stamp is included with the comment.
    list_emds : list
        A list of the EMD groups in the file.

    Examples
    --------
    Open an Berkeley EMD file using the *advanced* interface. See the emdReader function for a more convenient
    way to load the data. We show how to list the available data sets, load a 3D data set and plot the first image.
    >> from ncempy.io import emd
    >> with emd.fileEMD('filename.emd') as emd1:
    >>     [print(dataGroup.name) for dataGroup in emd1.list_emds] # print all available EMD datasets
    >>     data1, dims1 = emd1.get_emdgroup(0) # load the first full data array and dimension information
    """

    def __init__(self, filename, readonly=True):
        """Init opening/creating the file.

        Parameters
        ----------
        filename : str or pathlib.Path or file object
            The EMD file to open.
        readonly : bool, default True
            Set to False to allow writing to the file.

        """

        # necessary declarations in case something goes bad
        self.file_hdl = None

        # convenience handles to access the data in the emd file, everything can as well be accessed using the file_hdl
        self.version = None
        self.data = None
        self.microscope = None
        self.sample = None
        self.user = None
        self.comments = None
        self.list_emds = []  # list of HDF5 groups with emd_data_type type 1

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
            self.file_path = filename
            self.file_name = self.file_path.name

        # try opening the file
        if readonly:
            try:
                self.file_hdl = h5py.File(filename, 'r')
            except:
                print('Error opening file for readonly: "{}"'.format(filename))
                raise
        else:
            try:
                self.file_hdl = h5py.File(filename, 'a')
            except:
                print('Error opening file for read/write: "{}"'.format(filename))
                raise

        # if we got a working file
        if self.file_hdl:
            # check version information
            if 'version_major' in self.file_hdl.attrs and 'version_minor' in self.file_hdl.attrs:
                # read version information
                self.version = (self.file_hdl.attrs['version_major'], self.file_hdl.attrs['version_minor'])
            else:
                # set version information
                if not readonly:
                    self.file_hdl.attrs['version_major'] = 0
                    self.file_hdl.attrs['version_minor'] = 2

            # check for data group
            if 'data' not in self.file_hdl:
                if not readonly:
                    self.data = self.file_hdl.create_group('data')
            else:
                self.data = self.file_hdl['data']

            # check for data group
            if 'microscope' not in self.file_hdl:
                if not readonly:
                    self.microscope = self.file_hdl.create_group('microscope')
            else:
                self.microscope = self.file_hdl['microscope']

            # check for data group
            if 'sample' not in self.file_hdl:
                if not readonly:
                    self.sample = self.file_hdl.create_group('sample')
            else:
                self.sample = self.file_hdl['sample']

            # check for data group
            if 'user' not in self.file_hdl:
                if not readonly:
                    self.user = self.file_hdl.create_group('user')
            else:
                self.user = self.file_hdl['user']

            # check for data group
            if 'comments' not in self.file_hdl:
                if not readonly:
                    self.comments = self.file_hdl.create_group('comments')
            else:
                self.comments = self.file_hdl['comments']

            # find emd_data_type groups in the file
            self.list_emds = self.find_emdgroups(self.file_hdl)
            
            if len(self.list_emds) == 0 and readonly is True:
                message = 'No Berkeley EMD data sets found in file.'
                raise NoEmdDataSets(message)
            
            # compare to implementation
            # Only raise an exception if no Berkeley EMD data sets are found and the version does not match
            if self.version and readonly is True:
                if (len(self.list_emds) == 0) and (self.version != (0, 2)):
                    message = 'No Berkeley EMD data sets found. You are reading a version {}.{} EMD file, \
                               this implementation assumes version 0.2'.format(self.version[0], self.version[1])
                    raise UnsupportedEmdVersion(message)
            
    def __del__(self):
        """Destructor for EMD file object.

        """
        # close the file
        # if(not self.file_hdl.closed):
        self.file_hdl.close()

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

    def find_emdgroups(self, parent):
        """Find all emd_data_type groups within the group parent and return a list of references to their HDF5 groups.

        Parameters
        ----------
            parent: h5py.Group
                Handle to the parent group.

        Returns
        -------
            : list
                A list of h5py.Group handles to children groups
                being emd_data_type groups.

        """

        emds0 = []

        # recursive function to run and retrieve groups with emd_group_type set to 1
        def proc_group(group, emds):
            # take a look at each item in the group
            for item in group:
                # check if group
                if group.get(item, getclass=True) == h5py.Group:
                    item = group.get(item)
                    # check if emd_group_type
                    if 'emd_group_type' in item.attrs and 'data' in item.keys():
                        if item.attrs['emd_group_type'] == 1:
                            emds.append(item)
                    # process subgroups
                    proc_group(item, emds)

        # run
        proc_group(parent, emds0)

        return emds0

    def get_emddims(self, group):
        """Get the emdtype dimensions saved in in group.

        Parameters
        ----------
            group: h5py.Group
                Reference to the emdtype HDF5 group.

        Returns
        -------
            : tuple
                List of dimension vectors plus labels and units.

        """
        # get the dims
        dims = []
        for ii in range(group['data'].ndim):
            dim = group['dim{}'.format(ii + 1)]
            # save them as (vector, name, units)

            if 'name' in dim.attrs:
                if isinstance(dim.attrs['name'], np.ndarray):
                    name = dim.attrs['name'][0]
                else:
                    name = dim.attrs['name']
            else:
                name = 'dim{}'.format(ii + 1)

            if 'units' in dim.attrs:
                if isinstance(dim.attrs['units'], np.ndarray):
                    units = dim.attrs['units'][0]
                else:
                    units = dim.attrs['units']
            else:
                units = 'pixels'

            # Handle bytes objects by decoding them to strings
            # If something goes wrong, the original attribute is left as-is
            try:
                if isinstance(name, bytes):
                    name = name.decode('utf-8')
                if isinstance(units, bytes):
                    units = units.decode('utf-8')
            except:
                pass

            dims.append((dim[:], name, units))

        dims = tuple(dims)
        return dims

    def get_emdgroup(self, group):
        """Get the emd data saved in the requested group.

        Note
        ____
        The memmap keyword has been removed. Please use get_memmap().

        Parameters
        ----------
            group: h5py._hl.group.Group or int
                Reference to the HDF5 group to load. If int is used then the item corresponding to self.list_emds
                is loaded

        Returns
        -------
            : tuple or None
                None or tuple containing:
                    : np.ndarray
                        The data of the emdtype group.
                    : list
                        List of [0] dimension vectors, [1] labels and [2] units.
        """

        # check input
        if not isinstance(group, h5py.Group):
            if isinstance(group, int):
                try:
                    group = self.list_emds[group]
                except IndexError:
                    print('group does not exist')
                    return
            else:
                raise TypeError('group needs to refer to a valid HDF5 group!')

        if 'emd_group_type' not in group.attrs:
            raise TypeError('group is not a emd_group_type group!')
        if not group.attrs['emd_group_type'] == 1:
            raise TypeError('group is not a emd_group_type group!')

        # retrieve data and dims
        try:
            data = group['data'][:]
            dims = self.get_emddims(group)

            return data, dims
        except:
            # if something goes wrong, return None
            print('Content of "{}" does not seem to be in emd specified shape'.format(group.name))

            return None

    def get_memmap(self, group):
        """ Opens a new fileEMD object and returns the requested EMD group. This prevents the file from being closed.
        The original HDF5 file is left open. The file pointer is lost, but it will be closed once the EMD group
        is deleted or removed. See also get_emdgroup().


        Parameters
        ----------
        group: h5py._hl.group.Group or int
                Reference to the HDF5 group to load. If int is used then the item corresponding to self.list_emds
                is loaded

        """
        f = h5py.File(self.file_hdl.filename, 'r')

        # check input
        if not isinstance(group, h5py.Group):
            if isinstance(group, int):
                try:
                    group = self.list_emds[group]
                except IndexError:
                    print('group does not exist')
                    return
            else:
                raise TypeError('group needs to refer to a valid HDF5 group!')
        data = f[group.name + '/data']
        return data, self.get_emddims(group)

    def write_dim(self, label, dim, parent):
        """Auxiliary function to write a dim dataset to parent.

        Input is not checked for sanity, so handle exceptions in call.

        Parameters
        ----------
            label: str
                Label for dataset, usually dim1, dim2, dimN.
            dim: tuple
                Tuple containing (data, name, units).
            parent: h5py.Group
                HDF5 handle to parent group.

        Returns
        -------
            : h5py.Group
                HDF5 dataset handle referencing this dim.

        """

        try:
            dset = parent.create_dataset(label, data=dim[0])
            dset.attrs['name'] = np.string_(dim[1])
            dset.attrs['units'] = np.string_(dim[2])
        except:
            raise RuntimeError('Error during writing dim dataset')

        return dset

    def put_emdgroup(self, label, data, dims, parent=None, overwrite=False, **kwargs):
        """Put an emdtype dataset into the EMD file.

        Parameters
        ----------
            label: str
                Label for the emdtype group containing the dataset.
            data: np.ndarray
                Numpy array containing the data.
            dims: tuple
                Tuple containing the necessary dims as ((vec, name, units), (vec, name, units), ...)
            parent: h5py.Group or None
                Parent for the emdtype group, if None it will be written to /data.
            overwrite: bool
                Set to force overwriting entry in EMD file.
            **kwargs: various
                Keyword arguments to be passed to h5py.create_dataset(), e.g. for compression.

        Returns
        -------
            : h5py.Group or None
                Group referencing this emdtype dataset or None if failed.
        """

        # check input
        if not isinstance(label, str):
            raise TypeError('label needs to be string!')

        if not isinstance(data, np.ndarray):
            raise TypeError('data needs to be a numpy.ndarray!')

        try:
            assert len(dims) == len(data.shape)
            for i in range(len(dims)):
                assert len(dims[i]) == 3
                assert dims[i][0].shape[0] == data.shape[i]
        except:
            raise TypeError('Something wrong with the provided dims')

        # write stuff to HDF5

        # create group
        try:
            if parent:
                if label in parent:
                    if overwrite:
                        print('overwriting "{}" in "{}"'.format(label, parent.name))
                        del parent[label]
                    else:
                        print('"{}" already exists in "{}"'.format(label, parent.name))
                        raise RuntimeError('"{}" already exists in "{}"'.format(label, parent.name))
                grp = parent.create_group(label)

            else:
                if label in self.data:
                    if overwrite:
                        print('overwriting "{}" in "{}"'.format(label, self.data.name))
                        del self.data[label]
                    else:
                        print('"{}" already exists in "{}"'.format(label, self.data.name))
                        raise RuntimeError('"{}" already exists in "{}"'.format(label, self.data.name))

                grp = self.data.create_group(label)

            # add attribute
            grp.attrs['emd_group_type'] = 1

            # create dataset
            _ = grp.create_dataset('data', data=data, **kwargs)

            # create dim datasets
            for i in range(len(dims)):
                self.write_dim('dim{}'.format(i + 1), dims[i], grp)

            # update emds list
            self.list_emds = self.find_emdgroups(self.file_hdl)

            return grp

        except:
            print('Something went wrong trying to write the dataset.')

            return None

    def put_comment(self, msg, timestamp=None):
        """Create a comment in the EMD file.

        If timestamp already exists, the msg is appended to existing comment.

        Parameters
        ----------
            msg: str
                String of the message to save.
            timestamp: str/None
                Timestamp used as the key, defaults to the current UTC time.

        """

        # check input
        if not isinstance(msg, str):
            raise TypeError('msg needs to be a string!')

        # create timestamp if missing
        if not timestamp:
            timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S (UTC)')
        else:
            # try to convert given timestamp to string
            try:
                timestamp = str(timestamp)
            except:
                raise

        # write comment
        if timestamp in self.comments.attrs:
            # append to existing
            self.comments.attrs[timestamp] += np.string_('\n' + msg)

        else:
            # create new entry
            self.comments.attrs[timestamp] = np.string_(msg)
        

def defaultDims(data, pixel_size=None, pixel_unit=None):
    """ A helper function that can generate a properly setup dim tuple
    with default values to allow quick writing of EMD files without
    the need to create these dim vectors manually.

    Parameters
    ----------
        data : ndarray
            The data that will be written to the EMD file. This is used to get the number of dims and their shape
        pixel_size : tuple, optional
            A tuple of pixel sizes. Must have same length as the number of dimensions of data.
        pixel_unit : tuple, optional
            A tuple of pixel units (i.e. nm), Must have the same length as the number of dimensions of data.
    Returns
    -------
        dims: list
            A properly formatted tuple of dim vectors used as input
            to emd.emdFile.put_emdgroup()
    """
    if isinstance(pixel_size, np.ndarray):
        raise Exception('pixel_size must be a tuple')

    num = data.ndim

    if not pixel_size:
        pixel_size = (1,) * num
    
    if not pixel_unit:
        pixel_unit = [f'unit{ii+1}' for ii in range(num)]
    
    if len(pixel_size) != data.ndim or len(pixel_unit) != data.ndim:
        raise ValueError('pixel_size, pixel_unit, and data dimensions must match')

    dims = []
    for ii in range(num):
        curDim = [np.linspace(0, data.shape[ii] - 1, data.shape[ii]) * pixel_size[ii],
                  'dim{}'.format(ii+1),
                  pixel_unit[ii]]
        dims.append(curDim)

    return dims


def emdReader(filename, dsetNum=0):
    """ A simple helper function to read in the data and metadata 
    in a structured format similar to the other ncempy readers.

    Note
    ----
        Note fully implemented yet. Work in progress.

    Parameters
    ----------
        filename : str or pathlib.Path
            The path to the file as a string.
        dsetNum : int
            The index of the data set to load.
            
    Returns
    -------
        : dict
            Data and metadata as a dictionary similar to other ncempy readers.

    Example
    -------
        Simply load all data and metadata from a data set in an EMD Velox file
            >> import ncempy.io as nio
            >> emd0 = nio.emd.emdReader('filename.emd', dsetNum = 0)

    """
    with fileEMD(filename, readonly=True) as emd0:
        d, dims = emd0.get_emdgroup(dsetNum)  # memmap must be false. File is closed
        out = {'data': d, 'filename': filename, 'pixelSize': []}

        # Get the group name
        name = emd0.list_emds[dsetNum].name.split('/')[-1]
        out['name'] = name

        for dim in dims:
            try:
                d = dim[0][1] - dim[0][0]
            except:
                d = 0
            out['pixelSize'].append(d)
        out['pixelUnit'] = [aa[2] for aa in dims]
        out['pixelName'] = [aa[1] for aa in dims]
        return out


def emdWriter(filename, data, pixel_size=None, pixel_unit=None, overwrite=False):
    """ Simple method to write data to a file formatted as an EMD v0.2. The only possible metadata to write is the pixel
    size for each dimension. Use the emd.fileEMD() class for more complex operations. The file must not already exist.

    Parameters
    ----------
    filename : str or pathlib.Path
        The path and file name to write the data to.
    data : ndarray
        The data as an ndarray to write to the file.
    pixel_size : tuple
        A tuple with the same length as the number of dims in data. len(pixel_size) == data.ndim
    pixel_unit : tuple, optional
            A tuple of pixel units (i.e. nm), Must have the same length as the number of dimensions of data.
    overwrite : boolean
        If file exists, overwrite it.
    """
    if isinstance(filename, str):
        filename = Path(filename)
    
    if pixel_size:
        try:
            assert len(pixel_size) == data.ndim
        except ValueError:
            raise ValueError('pixel_size length must match the number of dimensions of data.')
    
    if pixel_unit:
        try:
            assert len(pixel_unit) == data.ndim
        except ValueError:
            raise ValueError('pixel_unit length must match the number of dimensions of data.')
    
    # Delete the file if it exists and overwrite is true
    if filename.exists() and overwrite:
        filename.unlink()

    if not filename.exists():
        with fileEMD(filename, readonly=False) as emd0:
            # Setup the dims for the EMD file. Pixel size is set to 1 by default for each dimension
            dims0 = defaultDims(data, pixel_size=pixel_size, pixel_unit=pixel_unit)
            # Write the data to he emd file.
            emd0.put_emdgroup('converted', data, dims0)
    else:
        raise FileExistsError
