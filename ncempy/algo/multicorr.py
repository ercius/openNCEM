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

    def initial_correlation_image(self):
        '''Generate correlation image at initial resolution'''
        G12 = np.multiply(self.G1, np.conj(self.G2)) # is this the correct order that we want?
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
        self.xyShift = list(np.unravel_index(self.imageCorrIFT.argmax(), self.imageCorrIFT.shape, 'C'))

        if self.upsampleFactor == 1:
            imageSize = self.imageCorrIFT.shape
            self.xyShift[0] = ((self.xyShift[0] + imageSize[0]/2) % imageSize[0]) - imageSize[0]/2
            self.xyShift[1] = ((self.xyShift[1] + imageSize[1]/2) % imageSize[1]) - imageSize[1]/2
            # print(self.xyShift)
            self.G2shift = np.fft.fft2(np.roll(np.roll(np.fft.ifft2(self.G2), int(self.xyShift[0]), 0), int(self.xyShift[1]), 1))
        else:
            self.imageCorrLarge = self.upsampleFFT(self.imageCorr, 2)
            imageSizeLarge = self.imageCorrLarge.shape
            xySubShift2 = list(np.unravel_index(self.imageCorrLarge.argmax(), imageSizeLarge, 'C'))
            xySubShift2[0] = ((xySubShift2[0] + imageSizeLarge[0]/2) % imageSizeLarge[0]) - imageSizeLarge[0]/2
            xySubShift2[1] = ((xySubShift2[1] + imageSizeLarge[1]/2) % imageSizeLarge[1]) - imageSizeLarge[1]/2
            self.xyShift = [i/2 for i in xySubShift2]
            print(['self.xyShift'] + self.xyShift)

            if self.upsampleFactor > 2:
                # here is where we use DFT registration to make things much faster
                # we cut out and upsample a peak 1.5 by 1.5 px from or ouriginal correlation image.

                self.xyShift[0] = np.round(self.xyShift[0] * self.upsampleFactor) / self.upsampleFactor
                self.xyShift[1] = np.round(self.xyShift[1] * self.upsampleFactor) / self.upsampleFactor

                globalShift = np.fix(np.ceil(self.upsampleFactor * 1.5)/2)# this line might have an off by one error based. The associated matlab comment is "this will be used to center the output array at dftshift + 1"
                print('globalShift', globalShift)

                self.imageCorrUpsample = np.conj(self.dftUpsample(np.conj(self.imageCorr), self.upsampleFactor, globalShift - np.multiply(self.xyShift, self.upsampleFactor))) / (np.fix(imageSizeLarge[0]) * np.fix(imageSizeLarge[1]) * self.upsampleFactor ** 2)

                self.xySubShift = np.unravel_index(self.imageCorrUpsample.argmax(), self.imageCorrUpsample.shape, 'C')
                print('self.xySubShift', self.xySubShift)

                # add a subpixel shift via parabolic fitting
                try:
                    icc = np.real(imageCorrUpsample[self.xySubShift[0] - 1 : self.xySubShift[0] + 2, self.xySubShift[1] - 1 : self.xySubShift[1] + 2])
                    dx = (icc[2,1] - icc[0,1]) / (4 * icc[1,1] - 2 * icc[2,1] - 2 * icc[0,1])
                    dy = (icc[1,2] - icc[1,0]) / (4 * icc[1,1] - 2 * icc[1,2] - 2 * icc[1,0])
                except:
                    dx, dy = 0, 0 # this is the case when the peak is near the edge and one of the above values does not exist

                self.xySubShift = self.xySubShift - globalShift;
                self.xyShift = self.xyShift + (self.xySubShift + np.array([dx, dy])) / self.upsampleFactor

    def imageShifter(self):
        imageSize = self.imageCorrIFT.shape
        qx = self.makeFourierCoords(imageSize[0], 1) # does this need to be a column vector
        if imageSize[1] == imageSize[0]:
            qy = qx
        else:
            qy = self.makeFourierCoords(imageSize[1], 1)

        self.G2shift = np.multiply(G2, np.multiply( np.exp(-2j * np.pi * qx * xyShift[0]),  np.exp(-2j * np.pi * qy * xyShift[1])))


    def upsampleFFT(self, imageInit, upsampleFactor):
        '''This does a fourier upsample of the imageInit. imageInit is the fourier transform of the correlation image. upsampleFactor is self-descriptive. The function returns the real space correlation image that has been fourier upsampled by 2.'''
        imageSize = imageInit.shape
        imageUpsample = np.zeros(tuple((i*upsampleFactor for i in imageSize))) + 0j
        imageUpsample[:imageSize[0], :imageSize[1]] = imageInit
        imageUpsample = np.roll(np.roll(imageUpsample, int(imageSize[0]/2), 0), int(imageSize[1]/2),1)
        imageUpsampleReal = np.real(np.fft.ifft2(imageUpsample))

        return imageUpsampleReal

    def makeFourierCoords(N, pSize):
        if N % 2 == 0:
            q = np.roll(np.arange(-N/2, N/2) / ((N-1) * pSize), -N/2, axis=1)
        else:
            q = np.roll(np.arange(-N/2 + 0.5, N/2 + 0.5) / ((N-1) * pSize), -N/2 + 0.5, axis=1)
        return q

    def dftUpsample(self, imageCorr, upsampleFactor, xyShift):
        '''this code is adapted from dftups which is a subfunction available from dftregistration on the matlab file exchange. At some point I will credit the people who did that code. It has been hardcoded with a 1.5 pixel'''
        imageSize = imageCorr.shape
        pixelRadius = 1.5
        numRow = np.ceil(pixelRadius * upsampleFactor)
        numCol = numRow
        # compute kernels and compute subpixel DFT by matrix multiplication. see here: http://www.sciencedirect.com/science/article/pii/S0045790612000778

        colKern = np.exp(
        (-1j * 2 * np.pi / (imageSize[1] * upsampleFactor))
        * (np.fft.ifftshift( (np.arange(imageSize[1])) )
        - np.floor(imageSize[1]/2))
        * (np.arange(numCol) - xyShift[1])[:, np.newaxis]
        ) # this is currently not a column, and it needs to be one for the matrix multiplication to work correctly. Something needs to be changed with the transpose portion, or the vector multiplication needs to be done differently.

        rowKern = np.exp(
        (-1j * 2 * np.pi / (imageSize[0] * upsampleFactor))
        * (np.arange(numRow) - xyShift[0])
        * (np.fft.ifftshift(np.arange(imageSize[0]))
        - np.floor(imageSize[0]/2))[:, np.newaxis]
        )
        print(['rowkern, imagecorr, colkern'] + [rowKern.shape, imageCorr.shape, colKern.shape])
        return np.real(np.dot(np.dot(rowKern.transpose(), imageCorr), colKern.transpose()))



    def multicorr(self):
        self.parse_input()
        self.initial_correlation_image()
        self.upsampled_correlation()
        # return self.xyShift, self.G2shift, self.imageCorrIFT
