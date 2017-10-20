import datetime
import inspect
import ncempy.io.dm
import os
import unittest


def get_measure(start=None, operation="Unk"):
    if start is None:
        return  datetime.datetime.now()
    else:
        end = datetime.datetime.now()
        delta = end-start
        print("PID({})-Operation ({}): {}".format(os.getpid(), operation,
                                                  delta))
        return delta
        
_curren_file = os.path.dirname(os.path.abspath(
                    inspect.getfile(inspect.currentframe())))
_images_folder="{}/resources".format(_curren_file)

def _get_image_route(file_name):
    
    return os.path.join(_images_folder, file_name)

class test_dm3(unittest.TestCase):
    '''
    Test the DM3 io module
    '''
    
    def _read_dm3_data(self, file_route, on_memory=False):
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
                                _get_image_route("Frame_115K_spot4_300V.dm3"))
        
        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    
    def test_read_dm3_on_memory(self):
        
        metadata, img = self._read_dm3_data(
                                _get_image_route("Frame_115K_spot4_300V.dm3"),
                                on_memory=True)
        
        
        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["metadata"]["pixelOrigin"], [0.0, 0.0])
        self.assertEqual(metadata["header"]['.ImageSourceList.1.ImageRef'], 1)
        
    def test_dm4_memory_vs(self):
        m0=get_measure()
        #file_route="/project/projectdirs/m2657/NCEM/TEAM_I/Jihan/FocalSeriesImages_tip3.dm4"
        metadata, img = self._read_dm3_data(
                                _get_image_route("FocalSeriesImages_tip3.dm4"))
        #metadata, img = self._read_dm3_data(file_route)
        delta0=get_measure(m0, "No buffer read")
        
        m1=get_measure()
        metadata, img = self._read_dm3_data(
                            _get_image_route("FocalSeriesImages_tip3_copy.dm4"),
                            on_memory=True)
        #metadata, img = self._read_dm3_data(file_route+".copy",
         #                   on_memory=True)
        delta1=get_measure(m1, "With buffer read")
        
        self.assertGreater(delta0, delta1)
        
        
