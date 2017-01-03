'''
Tests for the ser io module.
'''

import unittest
import ncempy.io.ser
import numpy as np
import os
import os.path

class test_ser(unittest.TestCase):
    '''
    Test the SER io module.
    '''

    def test_read_ser(self):
        '''
        Test the ser reading functionality.
        '''
        
        # wrong argument type
        with self.assertRaises(TypeError):
            fser = ncempy.io.ser.fileSER(42)

        # non existing file
        with self.assertRaises(IOError):
            fser = ncempy.io.ser.fileSER('')
            
        # wrong file
        with self.assertRaises(Exception):
            fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01.emi', verbose=True)

        # single 2D image file
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser')
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser', verbose=True)
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser', 'ncempy/test/resources/Pt_SAED_D910mm_single/im01.emi', verbose=True)
        
        # time series of 2D images
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01_1.ser')
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01_1.ser', verbose=True)
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01_1.ser','ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01.emi', verbose=True)
        
        # mapping of 2D images
        #fser = ncempy.io.ser.fileSER('ncempy/test/resources/CBED_mapping/map01_2.ser', verbose=True)
        
        # mapping of 1D datasets
        #fser = ncempy.io.ser.fileSER('ncempy/test/resources/CBED_mapping/map01_1.ser', verbose=True)
        
        # single 1D dataset
        # time series of 1D datasets
        # 2D mapping of 2D datasets
        # 2D mapping of 1D dataset
        ## not implemented yet
        
        
    def test_read_emi(self):
        '''
        Test the emi reading functionality.
        '''
        
        # ser file for this testing
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser')
        
        # wrong argument
        with self.assertRaises(TypeError):
            fser.readEMI(42)
        
        # non existing file
        with self.assertRaises(IOError):
            fser.readEMI('')
        
        # auxiliary function
        self.assertIsInstance(fser.parseEntryEMI('42'), int)
        self.assertIsInstance(fser.parseEntryEMI('42.42'), float)
        self.assertIsInstance(fser.parseEntryEMI('forty two'), np.string_)
        
        # read 
        fser.readEMI('ncempy/test/resources/Pt_SAED_D910mm_single/im01.emi')
        
        
    def test_read_dataset(self):
        '''
        Test functions for retrieving datasets.
        '''
    
        # ser file for this testing
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser')
        
        # wrong index
        with self.assertRaises(IndexError):
            fser.getDataset(-1)
            
        # wrong index type
        with self.assertRaises(TypeError):
            fser.getDataset('foo')
        
        # try dataset + meta
        dataset, meta = fser.getDataset(0, verbose=True)
        
        # try tag
        tag = fser.getTag(0, verbose=True)
        
    
    def test_badtagoffset(self):
        '''
        Bad TagOffsetArray found in some files.
        '''
        
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/ser_badtagoffset/01_Si110_5images_1.ser', 'ncempy/test/resources/ser_badtagoffset/01_Si110_5images.emi', verbose=True)
        
        for i in range(fser.head['ValidNumberElements']):
            tag = fser.getTag(i)
            
            assert(tag['TagTypeID'] == 0)
            assert(tag['Time'] == 0)
            assert(np.isnan(tag['PositionX']))
            assert(np.isnan(tag['PositionY']))
        
        if os.path.isfile('ncempy/test/resources/output/01_Si110_5images.emd'):
            os.remove('ncempy/test/resources/output/01_Si110_5images.emd')
        fser.writeEMD('ncempy/test/resources/output/01_Si110_5images.emd')
        
        #import pdb; pdb.set_trace()
        
        
        
    def test_write_emd(self):
        '''
        Test the emd writing functionality.
        '''
        
        # single 2D image file
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser')
        if os.path.isfile('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd'):
            os.remove('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd')
        fser.writeEMD('ncempy/test/resources/output/Pt_SAED_D910mm_single_woemi.emd')
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Pt_SAED_D910mm_single/im01_1.ser', 'ncempy/test/resources/Pt_SAED_D910mm_single/im01.emi')
        if os.path.isfile('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd'):
            os.remove('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd')
        fser.writeEMD('ncempy/test/resources/output/Pt_SAED_D910mm_single.emd')
        with self.assertRaises(IOError):
            fser.writeEMD('')
        
        # time series of 2D images
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01_1.ser', 'ncempy/test/resources/Au_SAED_D910mm_20x_at_800/pos01.emi')
        if os.path.isfile('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd'):
            os.remove('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd')
        fser.writeEMD('ncempy/test/resources/output/Au_SAED_D910mm_20x_at_800.emd')
        
        ## exceeds memory
        # mapping of 2D images
        #fser = ncempy.io.ser.fileSER('ncempy/test/resources/CBED_mapping/map01_2.ser')#, 'ncempy/test/resources/CBED_mapping/map01.emi')
        #if os.path.isfile('ncempy/test/resources/output/CBED_mapping_images.emd'):
        #    os.remove('ncempy/test/resources/output/CBED_mapping_images.emd')
        #fser.writeEMD('ncempy/test/resources/output/CBED_mapping_images.emd')
     
        
        # mapping of 1D datasets
        #fser = ncempy.io.ser.fileSER('ncempy/test/resources/CBED_mapping/map01_1.ser', 'ncempy/test/resources/CBED_mapping/map01.emi')
        #if os.path.isfile('ncempy/test/resources/output/CBED_mapping_edx.emd'):
        #    os.remove('ncempy/test/resources/output/CBED_mapping_edx.emd')
        #fser.writeEMD('ncempy/test/resources/output/CBED_mapping_edx.emd')
        
        
        # large time series of 2D images
        fser = ncempy.io.ser.fileSER('ncempy/test/resources/Au_SAED_D910mm_100x_at_RT/step_off_1.ser','ncempy/test/resources/Au_SAED_D910mm_100x_at_RT/step_off.emi', verbose=True)
        ##fser.head['ValidNumberElements'] = 20
        if os.path.isfile('ncempy/test/resources/output/Au_SAED_D910mm_100x_at_RT.emd'):
            os.remove('ncempy/test/resources/output/Au_SAED_D910mm_100x_at_RT.emd')
        fser.writeEMD('ncempy/test/resources/output/Au_SAED_D910mm_100x_at_RT.emd')
        
        # single 1D dataset
        # time series of 1D datasets
        # 2D mapping of 2D datasets
        # 2D mapping of 1D dataset
        ## not implemented yet


# to test with unittest runner
if __name__ == '__main__':
    unittest.main()
