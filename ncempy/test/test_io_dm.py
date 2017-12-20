'''
Tests for the basic functionalities of the dm io module.
'''

import datetime
import inspect
import ncempy.io.dm
import os
import unittest
from matplotlib.image import imread


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
        
        img = img3D
        if dimensions == 3:
            img = img[:,:,int(img.shape[2]/2),]
        elif dimensions == 4:
            img = img[:,:,
                  int(img.shape[2]/2),
                  int(img.shape[3]/2)]
        elif dimensions>4:
            raise ValueError("Images with more than four dimensions not"
                            " supported.")
        metadata = dict(dimensions=dimensions,
                              header=f.allTags,
                              metadata={x:ds[x]
                                        for x in ds.keys()
                                        if x!="data"})
        del f
        return metadata, img
    
    def test_read_dm3(self):
        
        metadata, img = self._read_dm3_data(self._get_image_route(
                            "dmTest_3D_float32_nonSquare_diffPixelSize.dm3"))

        self.assertEqual(metadata["dimensions"], 3)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    def test_read_dm3_one_dimension(self):
        
        metadata, img = self._read_dm3_data(self._get_image_route(
                            "06_lowLoss.dm3"))

        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 200.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)

        
    
    def test_read_dm3_on_memory(self):
        
        metadata, img = self._read_dm3_data(self._get_image_route(
                            "dmTest_3D_float32_nonSquare_diffPixelSize.dm3"),
                            on_memory=True)
        
        
        self.assertEqual(metadata["dimensions"], 3)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    def test_dm4_memory_vs_file_performance(self):
        """ Even with a local HD, memory read should be x10 faster."""
        m0=get_measure()
        for i in range(10):
            metadata, img = self._read_dm3_data(
                             self._get_image_route(
                                 "Si-SiGe-test-01-31x12x448x480.dm4"))
            delta0=get_measure(m0)
            
            m1=get_measure()
            metadata, img = self._read_dm3_data(
                                self._get_image_route(
                                    "Si-SiGe-test-01-31x12x448x480_copy.dm4"),
                                on_memory=True)
            delta1=get_measure(m1)
            
            self.assertGreater(delta0, delta1)
        
    def test_extract_on_memory(self):
        from matplotlib.image import imsave
        from matplotlib import cm
        f = ncempy.io.dm.fileDM(
            self._get_image_route("Si-SiGe-test-01-31x12x448x480.dm4"),
            on_memory=False)
        
        f.parseHeader()
        ds = f.getDataset(0)
        img3D_no_on_mem = ds['data']
        
        del f
        
        f = ncempy.io.dm.fileDM(
            self._get_image_route("Si-SiGe-test-01-31x12x448x480.dm4"),
            on_memory=True)
        
        f.parseHeader()
        ds = f.getDataset(0)
        img3D_on_mem = ds['data']
      
        del f
        self.assertTrue((img3D_no_on_mem==img3D_on_mem).all())
        
    def text_compare_png(self):
        
        for file_name in ["dmTest_3D_float32_nonSquare_diffPixelSize.dm3",
                          "dmTest_uint16_nonSquare_diffPixelSize.dm3",
                          "dmTest_int32_nonSquare_diffPixelSize.dm4"]:
            metadata, img = self._read_dm3_data(self._get_image_route(
                                "{}".format(file_name)),
                                on_memory=True)
            png = imread(self._get_image_route(
                            "{}.png".format(file_name)))
            self.assertEaual(img, png)
        
        