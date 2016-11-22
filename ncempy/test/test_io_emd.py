'''
Tests for the emd io module.
'''

import unittest
import os
import os.path
import numpy as np

import ncempy.io.emd

class test_emd(unittest.TestCase):
    '''
    Test the EMD io module
    '''
    
    def test_init_emd(self):
        
        # wrong argument type
        with self.assertRaises(TypeError):
            femd = ncempy.io.emd.fileEMD(42)
        
        # non existing file in readonly
        with self.assertRaises(IOError):
            femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/doesnotexist.emd', readonly=True)
            
        # impossible file for read/write
        with self.assertRaises(IOError):
            femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/output/output.emd')
            
        # open existing file
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/Au_SAED_D910mm_20x_at_800.emd')
        # open existing file readonly
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/Au_SAED_D910mm_20x_at_800.emd', readonly=True)
        # open nonexisting file for writing
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/test.emd')


    def test_data_manipulation(self):
    
        # open testfiles
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/Au_SAED_D910mm_20x_at_800.emd', readonly=True)
        if os.path.isfile('ncempy/test/resources/output/copy.emd'):
            os.remove('ncempy/test/resources/output/copy.emd')
        femd2 = ncempy.io.emd.fileEMD('ncempy/test/resources/output/copy.emd')
        
        # test error handling in get_emdgroup
        with self.assertRaises(TypeError):
            data, dims = femd.get_emdgroup(femd.file_hdl)
        
        # working
        data, dims = femd.get_emdgroup(femd.list_emds[0])
        
        # test error handling in put_emdgroup
        with self.assertRaises(TypeError):
            femd2.put_emdgroup(42, data,  dims)
            
        with self.assertRaises(TypeError):
            femd2.put_emdgroup('Au_SAED_D910mm_20x_at_800', 42, dims)
        
        with self.assertRaises(TypeError):
            femd2.put_emdgroup('Au_SAED_D910mm_20x_at_800', data, dims[0:1])
        
        # working
        self.assertIsNotNone(femd2.put_emdgroup('Au_SAED_D910mm_20x_at_800', data, dims))
        
        # try to write a readonly file
        femd2 = ncempy.io.emd.fileEMD('ncempy/test/resources/output/copy.emd', readonly=True)
        self.assertIsNone(femd2.put_emdgroup('Au_SAED_D910mm_20x_at_800', data, dims))
        
        del femd, femd2
        
        # write a testfile
        if os.path.isfile('ncempy/test/resources/output/test2.emd'):
            os.remove('ncempy/test/resources/output/test2.emd')
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/test2.emd')
        data = np.random.rand(512,512,100)
        dims = ( (np.array(range(512)), 'x', '[px]'),
                 (np.array(range(512)), 'y', '[px]'),
                 (np.linspace(0.0, 3.14, num=100), 'angle','[rad]') )
        self.assertIsNotNone(femd.put_emdgroup('dataset_1', data, dims))
        
        dim = (np.array(range(100)), 'number', '[]')
        femd.write_dim('dim3_number', dim, femd.list_emds[0])
        
        femd.put_comment('file created, filled with random numbers')
        
        # try to overwrite
        self.assertIsNone(femd.put_emdgroup('dataset_1', data, dims))
        
        
    def test_comments(self):
    
        # create a file for comments
        if os.path.isfile('ncempy/test/resources/output/comments.emd'):
            os.remove('ncempy/test/resources/output/comments.emd')
        femd = ncempy.io.emd.fileEMD('ncempy/test/resources/output/comments.emd')
        
        with self.assertRaises(TypeError):
            femd.put_comment(42)
        
        femd.put_comment('file created')
        femd.put_comment('how long does it take until the second comment is issued?')
        
        femd.put_comment('something happened', 'today')
        femd.put_comment('even more happened', 'today')

# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
