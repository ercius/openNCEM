'''
Tests for the basic functionalities of the dm io module.
'''

import datetime
import inspect
import ncempy.io.dm
import os
import unittest


def get_measure(start=None, operation=None):
    '''Auxiliary function to measure time spans.
    
    Args:
        start: if set to a datetime object, function returns a timespan
          object measurring the difference between start and now().
        operation: if set to a string, (and start is a datetime object),
          the function prints a message indicating the process PID, the
          string and the time_sona value.
    
    Returns: a datetime object (now()) or a timespan object.    
    '''
    if start is None:
        return  datetime.datetime.now()
    else:
        end = datetime.datetime.now()
        delta = end-start
        if operation:
            print("PID({})-Operation ({}): {}".format(os.getpid(), operation,
                                                  delta))
        return delta

         

class test_dm3(unittest.TestCase):
    '''
    Test the DM3 io module
    '''
    
    ''' Auxiliary attributes and method to facilitate the location of the test
    images'''
    _curren_file = os.path.dirname(os.path.abspath(
                    inspect.getfile(inspect.currentframe())))
    _images_folder="{}/resources".format(_curren_file)

    def _get_image_route(self, file_name):
        return os.path.join(self._images_folder, file_name)
    
    def _read_dm3_data(self, file_route, on_memory=False):
        '''Creates a DMobject and reads its data metadata
        
        Args:
             file_route: full route where the file is placed.
             on_memory: if True, the dm file will be opened in on memory mode.
            
        
        '''
        
        f = ncempy.io.dm.fileDM(file_route, on_memory=on_memory)
        if on_memory:
            self.assertTrue(f._on_memory)
        f.parseHeader()
        ds = f.getDataset(0)
        img3D = ds['data']
        dimensions=len(img3D.shape);
        if dimensions==2:
            img = img3D
        else:
            img = img3D[int(img3D.shape[0]/2),:,:]
        metadata = dict(dimensions=dimensions,
                              header=f.allTags,
                              metadata={x:ds[x]
                                        for x in ds.keys()
                                        if x!="data"})
        del f
        return metadata, img
    
    def test_read_dm3(self):
        
        metadata, img = self._read_dm3_data(
                            self._get_image_route("Frame_115K_spot4_300V.dm3"))
        
        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    
    def test_read_dm3_on_memory(self):
        
        metadata, img = self._read_dm3_data(
                            self._get_image_route("Frame_115K_spot4_300V.dm3"),
                            on_memory=True)
        
        
        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    def test_dm4_memory_vs_file(self):
        """ Even with a local HD, memory read should be x10 faster."""
        m0=get_measure()
        for i in range(10):
            metadata, img = self._read_dm3_data(
                             self._get_image_route(
                                 "FocalSeriesImages_tip3.dm4"))
            delta0=get_measure(m0)
            
            m1=get_measure()
            metadata, img = self._read_dm3_data(
                                self._get_image_route(
                                    "FocalSeriesImages_tip3_copy.dm4"),
                                on_memory=True)
            delta1=get_measure(m1)
            
            self.assertGreater(delta0, delta1)
        
        
