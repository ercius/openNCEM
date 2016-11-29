'''
This module provides an interface to the EMD file format.

See https://emdatasets.com/ for more details.
'''

import numpy as np
import h5py
import datetime


class fileEMD:
    '''Class to represent EMD files. 
    
    Implemented for spec 0.2 using the recommended layout for metadata.
    
    Meant to provide convenience functions for commonly occuring tasks. This means that you will still want to acces fileEMD.file_hdl to manipulate the HDF5 file for not so commonly occuring tasks.
    
    Parameters:
        filename (str):    Name of the EMD file.
        readonly (bool):    Set to open in read only mode.
            
    '''
    
    def __init__(self, filename, readonly=False):
        '''Init opening/creating the file.
        
        '''
        
        ## necessary declarations in case something goes bad
        self.file_hdl = None
        
        # convenience handles to access the data in the emd file, everything can as well be accessed using the file_hdl
        self.version = None
        self.data = None
        self.microscope = None
        self.sample = None
        self.user = None
        self.comments = None
        self.list_emds = []             # list of HDF5 groups with emd_data_type type
        
        # check for string
        if not isinstance(filename, str):
            raise TypeError('Filename is supposed to be a string!')

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
                # compare to implementation
                if not self.version == (0,2):
                    print('WARNING: You are reading a version {}.{} EMD file, this implementation assumes version 0.2!'.format(self.version[0], self.version[1]))
            else:
                # set version information
                if not readonly:
                    self.file_hdl.attrs['version_major'] = 0
                    self.file_hdl.attrs['version_minor'] = 2
                
            # check for data group
            if not 'data' in self.file_hdl:
                if not readonly:
                    self.data = self.file_hdl.create_group('data')
            else:
                self.data = self.file_hdl['data']
                
            # check for data group
            if not 'microscope' in self.file_hdl:
                if not readonly:
                    self.microscope = self.file_hdl.create_group('microscope')
            else:
                self.microscope = self.file_hdl['microscope']
                
            # check for data group
            if not 'sample' in self.file_hdl:
                if not readonly:
                    self.sample = self.file_hdl.create_group('sample')
            else:
                self.sample = self.file_hdl['sample']
                
            # check for data group
            if not 'user' in self.file_hdl:
                if not readonly:
                    self.user = self.file_hdl.create_group('user')
            else:
                self.user = self.file_hdl['user']
                
            # check for data group
            if not 'comments' in self.file_hdl:
                if not readonly:
                    self.comments = self.file_hdl.create_group('comments')
            else:
                self.comments = self.file_hdl['comments']
                
            # find emd_data_type groups in the file
            self.list_emds = self.find_emdgroups(self.file_hdl)
            

    def __del__(self):
        '''Destructor for EMD file object.
        
        '''
        
        # close the file
        if(self.file_hdl):
            self.file_hdl.close()


    def find_emdgroups(self, parent):
        '''Find all emd_data_type groups within the group parent and return a list of references to their HDF5 groups.
        
        Parameters:
            parent (h5py._hl.group.Group):  Handle to the parent group.
            
        Returns:
            (list):    A list of h5py._hl.group.Group handles to children groups being emd_data_type groups.
            
        '''
        
        emds = []
        
        # recursive function to run and retrieve groups with emd_group_type set to 1
        def proc_group(group, emds):
            # take a look at each item in the group
            for item in group:
                # check if group
                if group.get(item, getclass=True) == h5py._hl.group.Group:
                    item = group.get(item)
                    # check if emd_group_type
                    if 'emd_group_type' in item.attrs:
                        if item.attrs['emd_group_type'] == 1:
                            emds.append(item)
                    # process subgroups
                    proc_group(item, emds)
        
        # run
        proc_group(parent, emds)
        
        return emds


    def get_emdgroup(self, group):
        '''Get the emdtype data saved in in group.
        
        Parameters:
            group (h5py._hl.group.Group):   Reference to the emdtype HDF5 group.
        
        Returns:
            (tuple/None):    None or tuple containing:
            
                np.ndarray:    The data of the emdtype group.
                
                list:    List of dimension vectors plus labels and units.
                
        '''
        
        # check input
        if not isinstance(group, h5py._hl.group.Group):
            raise TypeError('group needs to refer to a valid HDF5 group!')
            
        if not 'emd_group_type' in group.attrs:
            raise TypeError('group is not a emd_group_type group!')
        if not group.attrs['emd_group_type'] == 1:
            raise TypeError('group is not a emd_group_type group!')

        # retrieve data
        try:
            # get the data
            data = group['data'][:]
            
            # get the dims
            dims = []
            for i in range(len(data.shape)):
                dim = group['dim{}'.format(i+1)]
                # save them as (vector, name, units)
                
                if isinstance(dim.attrs['name'], np.ndarray):
                    name = dim.attrs['name'][0]
                else:
                    name = dim.attrs['name']
                
                if isinstance(dim.attrs['units'], np.ndarray):
                    units = dim.attrs['units'][0]
                else:
                    units = dim.attrs['units']
                    
                dims.append( (dim[:], name.decode('utf-8'), units.decode('utf-8')) )
            
            dims = tuple(dims)
            
            return data, dims
            
        except:
            # if something goes wrong, return None
            print('Content of "{}" does not seem to be in emd specified shape'.format(group.name))
            
            return None
            

    def write_dim(self, label, dim, parent):
        '''Auxiliary function to write a dim dataset to parent.
        
        Input is not checked for sanity, so handle exceptions in call.
        
        Parameters:
            label (str):    Label for dataset, usually dim1, dim2, dimN.
            dim (tuple):    Tuple containing (data, name, units).
            parent (h5py._hl.group.Group):    HDF5 handle to parent group.
        
        Returns:
            (h5py._hl.group.Group):    HDF5 dataset handle referencing this dim.
            
        '''
        
        try:
            dset = parent.create_dataset(label, data=dim[0])
            dset.attrs['name'] = np.string_(dim[1])
            dset.attrs['units'] = np.string_(dim[2])
        except:
            raise RuntimeError('Error during writing dim dataset')
        
        return dset
        
        
    def put_emdgroup(self, label, data, dims, parent=None, overwrite=False):
        '''Put an emdtype dataset into the EMD file.
        
        Parameters:
            label (str):    Label for the emdtype group containing the dataset.
            data (np.ndarray):    Numpy array containing the data.
            dims (tuple):    Tuple containing the necessary dims as ((vec, name, units), (vec, name, units), ...)
            parent (h5py._hl.group.Group/None):    Parent for the emdtype group, if None it will be written to /data.
            overwrite (bool):    Set to force overwriting entry in EMD file.
        
        Returns:
            (h5py._hl.group.Group/None):    Group referencing this emdtype dataset or None if failed.
            
        '''
        
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
            dset = grp.create_dataset('data', data=data)
             
            # create dim datasets
            for i in range(len(dims)):
                self.write_dim('dim{}'.format(i+1), dims[i], grp)
                    
            # update emds list
            self.list_emds = self.find_emdgroups(self.file_hdl)
                    
            return grp
           
        except:
            print('Something went wrong trying to write the dataset.')
                
            return None


    def put_comment(self, msg, timestamp=None):
        '''Create a comment in the EMD file.
        
        If timestamp already exists, the msg is appended to existing comment.
        
        Parameters:
            msg (str):    String of the message to save.
            timestamp (str/None):    Timestamp used as the key, defaults to the current UTC time.
        
        '''
        
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
            self.comments.attrs[timestamp] += np.string_('\n'+msg)
        
        else:
            # create new entry
            self.comments.attrs[timestamp] = np.string_(msg)
        

