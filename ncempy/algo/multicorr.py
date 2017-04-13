'''
Module to correlate two images, functionally written.
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

    Returns:
        xyShift - the shift between G1 and G2 in pixels

    Defaults -
        method = 'phase';
        peakSearchMode = false;  By default, we return only global best shift value.
        upsampleFactor = 1;
    '''
    method, upsampleFactor = parse_input(G1, G2, method, upsampleFactor)
    imageCorr = initial_correlation_image(G1, G2, method, upsampleFactor)
    xyShift = upsampled_correlation(imageCorr, upsampleFactor)
    print(method)
    print('------------')
    return xyShift

def parse_input(G1, G2, method = 'cross', upsampleFactor = 1):
    '''This function parses the inputs to make sure they're correct. Will raise errors if G1 or G2 are not ndarrays.

    Inputs:
        G1 - Fourier transform of image we are aligning to. Image reference
        G2 - Fourier transform of image that is being aligned. Image being registered
        method - correlation method ('phase', 'cross', 'hybrid'). Default 'cross'
        upsampleFactor - scalar integer specifying 1/subpixel_precision of fit. Default = 1.

    Returns:
        method - string ('phase', 'cross', 'hybrid')
        upsampleFactor - scalar integer
    '''
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
    '''
    Generate correlation image at initial resolution using the method specified.

    Inputs:
        G1 - Fourier transform of image we are aligning to. Image reference
        G2 - Fourier transform of image that is being aligned. Image being registered
        method - correlation method ('phase', 'cross', 'hybrid'). Default 'cross'
        upsampleFactor - scalar integer specifying 1/subpixel_precision of fit. Default = 1.

    Returns:
        imageCorr - an ndarray correlation image that has not yet been inverse Fourier transformed.
    '''
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
    '''
    Upsamples the correlation image by a set integer factor upsampleFactor. If upsampleFactor == 2, then it is naively Fourier upsampled. If the upsampleFactoris higher than 2, then it uses dftUpsample, which is a more efficient way to Fourier upsample the image.

    Inputs:
        imageCorr - Fourier transformed correlation image returned by initial_correlation_image. Is an ndarray.
        upsampleFactor - scalar integer of how much upsampling should be performed.

    Returns:
        xyShift - two element list with shift in x and y of G2 with respect to G1.
    '''

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
        print(['xySubShift2 '] + xySubShift2)
        xySubShift2[0] = ((xySubShift2[0] + imageSizeLarge[0]/2) % imageSizeLarge[0]) - imageSizeLarge[0]/2
        xySubShift2[1] = ((xySubShift2[1] + imageSizeLarge[1]/2) % imageSizeLarge[1]) - imageSizeLarge[1]/2
        xyShift = [i/2 for i in xySubShift2] #signs have to flip, or mod wrong?
        # print(xySubShift2)
        print(['xyShiftln127'] + xyShift)

        if upsampleFactor > 2:
            # here is where we use DFT registration to make things much faster
            # we cut out and upsample a peak 1.5 by 1.5 px from our original correlation image.

            xyShift[0] = np.round(xyShift[0] * upsampleFactor) / upsampleFactor
            xyShift[1] = np.round(xyShift[1] * upsampleFactor) / upsampleFactor

            globalShift = np.fix(np.ceil(upsampleFactor * 1.5)/2)# this line might have an off by one error based. The associated matlab comment is "this will be used to center the output array at dftshift + 1"
            print('globalShift', globalShift, 'upsampleFactor', upsampleFactor, 'xyShift', xyShift)

            imageCorrUpsample = np.conj(dftUpsample(np.conj(imageCorr), upsampleFactor, globalShift - np.multiply(xyShift, upsampleFactor))) / (np.fix(imageSizeLarge[0]) * np.fix(imageSizeLarge[1]) * upsampleFactor ** 2)

            xySubShift = np.unravel_index(imageCorrUpsample.argmax(), imageCorrUpsample.shape, 'C')
            # xySubShift = np.add(list(xySubShift), [1, 1])
            print('xySubShift', xySubShift)

            # add a subpixel shift via parabolic fitting
            try:
                icc = np.real(imageCorrUpsample[xySubShift[0] - 1 : xySubShift[0] + 2, xySubShift[1] - 1 : xySubShift[1] + 2])
                dx = (icc[2,1] - icc[0,1]) / (4 * icc[1,1] - 2 * icc[2,1] - 2 * icc[0,1])
                dy = (icc[1,2] - icc[1,0]) / (4 * icc[1,1] - 2 * icc[1,2] - 2 * icc[1,0])
            except:
                dx, dy = 0, 0 # this is the case when the peak is near the edge and one of the above values does not exist
            print('dxdy', dx, dy)
            print('xyShift', xyShift)
            xySubShift = xySubShift - globalShift;
            print('xysubShift2', xySubShift)
            xyShift = xyShift + (xySubShift + np.array([dx, dy])) / upsampleFactor
            print('xyShift2', xyShift)

    return xyShift

def upsampleFFT(imageInit, upsampleFactor):
    '''
    This does a Fourier upsample of the imageInit. imageInit is the Fourier transform of the correlation image. upsampleFactor is self-descriptive. The function returns the real space correlation image that has been Fourier upsampled by 2. It is written generally such that upsampleFactor can be greater than 2, but that should never happend/it has not been tested.

    The way it works is that it embeds imageInit in a larger array of zeros, then does the inverse Fourier transform to return the Fourier upsampled image in real space.

    Inputs:
        imageInit - ndarray of the image to be Fourier upsampled. This should be in the Fourier domain.
        upsampleFactor - integer scalar, almost always 2.

    Returns:
        imageUpsampleReal - the inverse Fourier transform of imageInit upsampled by the upsampleFactor. Is an ndarray.
    '''
    imageSize = imageInit.shape
    imageUpsample = np.zeros(tuple((i*upsampleFactor for i in imageSize))) + 0j
    imageUpsample[:imageSize[0], :imageSize[1]] = imageInit
    # plt.figure(1)
    # plt.imshow(np.real(imageUpsample))
    # plt.show(block = True)
    imageUpsample = np.roll(np.roll(imageUpsample, -int(imageSize[0]/2), 0), -int(imageSize[1]/2),1)
    imageUpsampleReal = np.real(np.fft.ifft2(imageUpsample))
    # plt.figure(1)
    # plt.imshow(imageUpsampleReal)
    # plt.show(block = True)
    return imageUpsampleReal

def dftUpsample(imageCorr, upsampleFactor, xyShift):
    '''
    This performs a matrix multiply DFT around a small neighboring region of the inital correlation peak. By using the matrix multiply DFT to do the Fourier upsampling, the efficiency is wildly improved. This is adapted from the subfuction dftups found in the dftregistration function on the Matlab File Exchange.

    https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

    The matrix multiplication DFT is from Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008). http://www.sciencedirect.com/science/article/pii/S0045790612000778

    Inputs:
        imageCorr - correlation image between two images in Fourier space. ndarray.
        upsampleFactor - scalar integer of how much to upsample.
        xyShift - single pixel shift between images previously computed. Used to center the matrix multiplication on the correlation peak. Is a two element list.
    Returns:
        imageUpsample - upsampled image from region around correlation peak. Is a ndarray (and I think the conjugate of the upsampled peak. Has to do wtih order of operations?)
    '''
    imageSize = imageCorr.shape
    pixelRadius = 1.5
    numRow = np.ceil(pixelRadius * upsampleFactor)
    numCol = numRow

    colKern = np.exp(
    (-1j * 2 * np.pi / (imageSize[1] * upsampleFactor))
    * (np.fft.ifftshift( (np.arange(imageSize[1])) )
    - np.floor(imageSize[1]/2))
    * (np.arange(numCol) - xyShift[1])[:, np.newaxis]
    ) # I think this can be written differently without the need for np.newaxis. This might require np.outer to compute the matrix itself instead of just using np.dot.

    rowKern = np.exp(
    (-1j * 2 * np.pi / (imageSize[0] * upsampleFactor))
    * (np.arange(numRow) - xyShift[0])
    * (np.fft.ifftshift(np.arange(imageSize[0]))
    - np.floor(imageSize[0]/2))[:, np.newaxis]
    ) # Comment from above applies.

    imageUpsample = np.real(np.dot(np.dot(rowKern.transpose(), imageCorr), colKern.transpose()))

    return imageUpsample

def imageShifter(G2, xyShift):
    '''
    This function multiplies G2 by a plane wave that has the real space effect of shifting ifft2(G2) by [x, y] pixels.

    Inputs:
        G2 - the Fourier transform of an image. ndarray.
        xyShift - a two element list

    Returns:
        G2shift - Fourier shifted G2. ndarray.
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
    '''
    This function creates Fourier coordinates such that (0,0) is in the center of the array.

    Inputs:
        N - the maximum coordinate in the original frame.
        pSize - the pixel size

    Returns:
        q - a single row array that has transformed 0:N to -N/2:N/2, such that the array sizes are the same.
    '''

    N = float(N)
    if N % 2 == 0:
        q = np.roll(np.arange(-N/2, N/2, dtype = 'float64') / (N * pSize), int(-N/2), axis=0)
    else:
        q = np.roll(np.arange((1-N)/2, (N+1)/2) / (N * pSize), int((1-N)/2), axis=0)
    return q
