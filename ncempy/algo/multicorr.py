'''
Module to correlate two images, functionally written.
'''

import numpy as np
import matplotlib
import cv2

def multicorr(G1, G2, method = 'cross', upsampleFactor = 1):
    '''
    Align a template to an image, possibly with multiple alignment peaks.
    Inputs:
        G1 - Fourier transform of image we are aligning to. Image reference
        G2 - Fourier transform of image that is being aligned. Image being registered
    Optional Inputs:
        'phase' or 'cross' or 'hybrid' - string that specifies correlation method used. (default = 'cross')
        upsampleFactor - scalar integer specifying 1/subpixel_precision of fit. (default = 1)

    Outputs:
        xyShift - the shift between G1 and G2 in pixels

    Defaults -
        method = 'phase';
        peakSearchMode = false;  By default, we return only global best shift value.
        upsampleFactor = 1;
    '''
    method, upsampleFactor = parse_input(G1, G2, method, upsampleFactor)
    imageCorr = initial_correlation_image(G1, G2, method, upsampleFactor)
    xyShift = upsampled_correlation(imageCorr, upsampleFactor)
    return xyShift

def parse_input(G1, G2, method = 'cross', upsampleFactor = 1):
    '''This function parses the inputs to make sure they're correct'''
    # Check to make sure both G1 and G2 are arrays
    if type(G1) is not np.ndarray:
        raise TypeError('G1 must be an ndarray')
    elif type(G2) is not np.ndarray:
        raise TypeError('G2 must be an ndarray')

    # Check to make sure method and upsample factor are the correct values
    if method not in ['phase', 'cross', 'hybrid']:
        print('Unknown method used, setting to cross')
        method = 'cross'

    if type(upsampleFactor) is not int and type(upsampleFactor) is not float:
        print('Upsample factor is not an integer or float, setting to 1')
        upsampleFactor = 1
    elif type(upsampleFactor) is not int:
        print('Upsample factor is not an integer, rounding down')
        upsampleFactor = int(upsampleFactor)
        if upsampleFactor < 1:
            print('Upsample factor is < 1, setting to 1')
            upsampleFactor = 1

    # Verify images are the same size.
    if G1.shape != G2.shape:
        raise TypeError('G1 and G2 are not the same size, G1 is {0} and G2 is {1}'.format(str(G1.shape), str(G2.shape)))

    return method, upsampleFactor

def initial_correlation_image(G1, G2, method = 'cross', upsampleFactor = 1):
    '''Generate correlation image at initial resolution'''
    G12 = np.multiply(G2, np.conj(G1)) # is this the correct order that we want?
    if method == 'phase':
        imageCorr = np.exp(1j * np.angle(G12))
    elif method == 'cross':
        imageCorr = G12
    elif method == 'hybrid':
        imageCorr = np.multiply(np.sqrt(np.absolute(G12)), np.exp(1j * np.angle(G12)))
    else:
        raise TypeError('{} method is not allowed'.format(str(method)))

    return imageCorr

def upsampled_correlation(imageCorr, upsampleFactor):
    imageCorrIFT = np.real(np.fft.ifft2(imageCorr))
    xyShift = list(np.unravel_index(imageCorrIFT.argmax(), imageCorrIFT.shape, 'C'))
    # print(['xyShift pre mod '] + xyShift)
    if upsampleFactor == 1:
        imageSize = imageCorrIFT.shape
        xyShift[0] = ((xyShift[0] + imageSize[0]/2) % imageSize[0]) - imageSize[0]/2
        xyShift[1] = ((xyShift[1] + imageSize[1]/2) % imageSize[1]) - imageSize[1]/2
        # print(['xyShift post mod '] + xyShift)
        #G2shift = np.fft.fft2(np.roll(np.roll(np.fft.ifft2(G2), int(xyShift[0]), 0), int(xyShift[1]), 1))
    else:
        imageCorrLarge = upsampleFFT(imageCorr, 2)
        imageSizeLarge = imageCorrLarge.shape
        xySubShift2 = list(np.unravel_index(imageCorrLarge.argmax(), imageSizeLarge, 'C'))
        # print(['xySubShift2 '] + xySubShift2)
        xySubShift2[0] = ((xySubShift2[0] + imageSizeLarge[0]/2) % imageSizeLarge[0]) - imageSizeLarge[0]/2
        xySubShift2[1] = ((xySubShift2[1] + imageSizeLarge[1]/2) % imageSizeLarge[1]) - imageSizeLarge[1]/2
        xyShift = [i/2 for i in xySubShift2] #signs have to flip, or mod wrong?
        # print(xySubShift2)
        # print(['xyShift'] + xyShift)

        if upsampleFactor > 2:
            # here is where we use DFT registration to make things much faster
            # we cut out and upsample a peak 1.5 by 1.5 px from our original correlation image.

            xyShift[0] = np.round(xyShift[0] * upsampleFactor) / upsampleFactor
            xyShift[1] = np.round(xyShift[1] * upsampleFactor) / upsampleFactor

            globalShift = np.fix(np.ceil(upsampleFactor * 1.5)/2)# this line might have an off by one error based. The associated matlab comment is "this will be used to center the output array at dftshift + 1"
            # print('globalShift', globalShift)

            imageCorrUpsample = np.conj(dftUpsample(np.conj(imageCorr), upsampleFactor, globalShift - np.multiply(xyShift, upsampleFactor))) / (np.fix(imageSizeLarge[0]) * np.fix(imageSizeLarge[1]) * upsampleFactor ** 2)

            xySubShift = np.unravel_index(imageCorrUpsample.argmax(), imageCorrUpsample.shape, 'C')
            # print('xySubShift', xySubShift)

            # add a subpixel shift via parabolic fitting
            try:
                icc = np.real(imageCorrUpsample[xySubShift[0] - 1 : xySubShift[0] + 2, xySubShift[1] - 1 : xySubShift[1] + 2])
                dx = (icc[2,1] - icc[0,1]) / (4 * icc[1,1] - 2 * icc[2,1] - 2 * icc[0,1])
                dy = (icc[1,2] - icc[1,0]) / (4 * icc[1,1] - 2 * icc[1,2] - 2 * icc[1,0])
            except:
                dx, dy = 0, 0 # this is the case when the peak is near the edge and one of the above values does not exist

            xySubShift = xySubShift - globalShift;
            xyShift = xyShift + (xySubShift + np.array([dx, dy])) / upsampleFactor

    return xyShift

def upsampleFFT(imageInit, upsampleFactor):
    '''This does a fourier upsample of the imageInit. imageInit is the fourier transform of the correlation image. upsampleFactor is self-descriptive. The function returns the real space correlation image that has been fourier upsampled by 2.'''
    imageSize = imageInit.shape
    imageUpsample = np.zeros(tuple((i*upsampleFactor for i in imageSize))) + 0j
    imageUpsample[:imageSize[0], :imageSize[1]] = imageInit
    imageUpsample = np.roll(np.roll(imageUpsample, int(imageSize[0]/2), 0), int(imageSize[1]/2),1)
    imageUpsampleReal = np.real(np.fft.ifft2(imageUpsample))

    return imageUpsampleReal

def dftUpsample(imageCorr, upsampleFactor, xyShift):
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
    # print(['rowkern, imagecorr, colkern'] + [rowKern.shape, imageCorr.shape, colKern.shape])
    return np.real(np.dot(np.dot(rowKern.transpose(), imageCorr), colKern.transpose()))

def imageShifter(G2, xyShift):
    '''
    This function fourier shifts G2 by [x, y] pixels.
    G2 is an image
    xyShift is a two element list
    '''
    imageSize = G2.shape
    qx = makeFourierCoords(imageSize[0], 1) # does this need to be a column vector
    if imageSize[1] == imageSize[0]:
        qy = qx
    else:
        qy = makeFourierCoords(imageSize[1], 1)

    G2shift = np.multiply(G2, np.outer( np.exp(-2j * np.pi * qx * xyShift[0]),  np.exp(-2j * np.pi * qy * xyShift[1])))

    return G2shift

def makeFourierCoords(N, pSize):
    N = float(N)
    if N % 2 == 0:
        q = np.roll(np.arange(-N/2, N/2, dtype = 'float64') / (N * pSize), int(-N/2), axis=0)
    else:
        q = np.roll(np.arange(-N/2 + 0.5, N/2 + 0.5) / ((N-1) * pSize), int(-N/2 + 0.5), axis=0)
    return q
