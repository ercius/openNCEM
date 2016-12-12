'''
Module to correlate two images.
'''

import numpy as np
import matplotlib
import cv2

class multicorr(object):
    '''
    Align a template to an image, possibly with multiple alignment peaks.
    Inputs:
        G1 - Fourier transform of image we are aligning to. Image reference
        G2 - Fourier transform of image that is being aligned. Image being registered
    Optional Inputs:
        'phase' or 'cross' or 'hybrid' - string that specifies correlation method used. (default = 'cross')
        upsampleFactor - scalar integer specifying 1/subpixel_precision of fit. (default = 1)
        peakList - [N x 3] array consisting of the (rough) peak locations to search for, where each row consists of [x0 y0 searchRadius]. (default - return only highest global peak)
        Note that (x0,y0) are the shift vectors, where (0,0) is unshifted.
    Outputs:
        xyShift - the shift between G1 and G2 in pixels
        G2shift - the Fourier transform of G2 such that it is shifted to most closely correlate with G1
        imageCorrIFT - the correlogram for G1 and G2, useful for debugging.

    Defaults -
        method = 'phase';
        peakSearchMode = false;  By default, we return only global best shift value.
        upsampleFactor = 1;
    '''
    def __init__(self, G1, G2, method='cross', upsampleFactor=1):
        self.G1 = G1
        self.G2 = G2
        self.method = method
        self.upsampleFactor = upsampleFactor

    def parse_input(self):
        '''This function parses the inputs to make sure they're correct'''
        # Check to make sure both G1 and G2 are arrays
        if type(self.G1) is not np.ndarray:
            raise TypeError('G1 must be an ndarray')
        elif type(self.G2) is not np.ndarray:
            raise TypeError('G2 must be an ndarray')

        # Check to make sure method and upsample factor are the correct values
        if self.method not in ['phase', 'cross', 'hybrid']:
            print('Unknown method used, setting to cross')
            self.method = 'cross'

        if type(self.upsampleFactor) is not int and type(self.upsampleFactor) is not float:
            print('Upsample factor is not an integer or float, setting to 1')
            self.upsampleFactor = 1
        elif type(self.upsampleFactor) is not int:
            print('Upsample factor is not an integer, rounding down')
            self.upsampleFactor = int(self.upsampleFactor)
            if self.upsampleFactor < 1:
                print('Upsample factor is < 1, setting to 1')
                self.upsampleFactor = 1

        # Verify images are the same size.
        if self.G1.shape != self.G2.shape:
            raise TypeError('G1 and G2 are not the same size, G1 is {0} and G2 is {1}'.format(str(self.G1.shape), str(self.G2.shape)))

    def initial_correlation(self):
        '''Generate correlation image at initial resolution'''
        G12 = np.multiply(self.G2, np.conj(self.G1)) # is this the correct order that we want? 
        if self.method == 'phase':
            imageCorr = np.exp(1j * np.angle(G12))
        elif self.method == 'cross':
            imageCorr = G12
        elif self.method == 'hybrid':
            imageCorr = np.multiply(np.sqrt(np.absolute(G12)), np.exp(1j * np.angle(G12)))
        else:
            raise TypeError('{} method is not allowed'.format(str(method)))
        self.imageCorr = imageCorr

    def upsampled_correlation(self):
        self.imageCorrIFT = np.real(np.fft.ifft2(self.imageCorr))
        self.xyShift = np.unravel_index(self.imageCorrIFT.argmax(), self.imageCorrIFT.shape, 'C')

    def multicorr(self):
        self.parse_input()
        self.initial_correlation()
        self.upsampled_correlation()
        # return self.xyShift, self.G2shift, self.imageCorrIFT
