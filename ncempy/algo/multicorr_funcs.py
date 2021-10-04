"""
Module to correlate two images, functionally written.

TODO
----
    - Cant use rfft2 currently. This gives one shift as 1/2 the value. How
      can this be improved to improve speed?
"""

import numpy as np


def initial_correlation_image(g1, g2, method='cross', verbose=False):
    """Generate correlation image at initial resolution using the method specified.

    Parameters
    ----------
        g1 : complex ndarray
            Fourier transform of reference image.
        g2 : complex ndarray
            Fourier transform of the image to register (the kernel).
        method : str, optional
            The correlation method to use. Must be 'phase' or 'cross' or 'hybrid' (default = 'cross')
        verbose : bool, default is False
            Print output.
    Returns
    -------
        imageCorr : ndarray complex
            Correlation array which has not yet been inverse Fourier transformed.
    """
    if verbose:
        print('Method is {}'.format(method))
    G12 = g2 * np.conj(g1)  # is this the correct order that we want?
    if method == 'phase':
        imageCorr = np.exp(1j * np.angle(G12))
    elif method == 'cross':
        imageCorr = G12
    elif method == 'hybrid':
        imageCorr = np.sqrt(np.abs(G12)) * np.exp(1j * np.angle(G12))
    else:
        raise TypeError('{} method is not allowed'.format(str(method)))

    return imageCorr


def upsampled_correlation(image_corr, upsample_factor, verbose=False):
    """Upsamples the correlation image by a set integer factor upsample_factor.
    If upsample_factor == 2, then it is naively Fourier upsampled.
    If the upsample_factor is higher than 2, then it uses dftUpsample, which is
    a more efficient way to Fourier upsample the image.

    Parameters
    ----------
        image_corr : ndarray complex
            Fourier transformed correlation image returned by initial_correlation_image.
        upsample_factor : int
            Upsampling factor.
        verbose : bool
            Provide output for debugging

    Returns
    -------
        xyShift : list
            Shift in x and y of G2 with respect to G1.
    """
    imageCorrIFT = np.real(np.fft.ifft2(image_corr))
    xyShift = list(np.unravel_index(imageCorrIFT.argmax(), imageCorrIFT.shape, 'C'))
    if verbose:
        print('xyShift initial = {}'.format(xyShift))
    if upsample_factor == 1:
        imageSize = imageCorrIFT.shape
        xyShift[0] = ((xyShift[0] + imageSize[0]/2) % imageSize[0]) - imageSize[0]/2
        xyShift[1] = ((xyShift[1] + imageSize[1]/2) % imageSize[1]) - imageSize[1]/2
    else:
        imageCorrLarge = upsampleFFT(image_corr, 2)
        imageSizeLarge = imageCorrLarge.shape
        xySubShift2 = list(np.unravel_index(imageCorrLarge.argmax(), imageSizeLarge, 'C'))
        if verbose:
            print('xySubShift2 = {}'.format(xySubShift2))
        xySubShift2[0] = ((xySubShift2[0] + imageSizeLarge[0]/2) % imageSizeLarge[0]) - imageSizeLarge[0]/2
        xySubShift2[1] = ((xySubShift2[1] + imageSizeLarge[1]/2) % imageSizeLarge[1]) - imageSizeLarge[1]/2
        xyShift = [i/2 for i in xySubShift2]  # signs have to flip, or mod wrong?
        if verbose:
            print('xyShift line 127 = {}'.format(xyShift))

        if upsample_factor > 2:
            # here is where we use DFT registration to make things much faster
            # we cut out and upsample a peak 1.5 by 1.5 px from our original correlation image.

            xyShift[0] = np.round(xyShift[0] * upsample_factor) / upsample_factor
            xyShift[1] = np.round(xyShift[1] * upsample_factor) / upsample_factor

            globalShift = np.fix(np.ceil(upsample_factor * 1.5) / 2)
            if verbose:
                print('globalShift = {}'.format(globalShift))
                print('xyShift = {}'.format(xyShift))

            imageCorrUpsample = np.conj(dftUpsample(np.conj(image_corr), upsample_factor,
                                                    globalShift - np.multiply(xyShift, upsample_factor))) / (np.fix(imageSizeLarge[0]) * np.fix(imageSizeLarge[1]) * upsample_factor ** 2)

            xySubShift = np.unravel_index(imageCorrUpsample.argmax(), imageCorrUpsample.shape, 'C')
            if verbose:
                print('xySubShift = {}'.format(xySubShift))

            # add a subpixel shift via parabolic fitting
            try:
                icc = np.real(imageCorrUpsample[xySubShift[0] - 1 : xySubShift[0] + 2,
                                                xySubShift[1] - 1 : xySubShift[1] + 2])
                dx = (icc[2, 1] - icc[0, 1]) / (4 * icc[1, 1] - 2 * icc[2, 1] - 2 * icc[0, 1])
                dy = (icc[1, 2] - icc[1, 0]) / (4 * icc[1, 1] - 2 * icc[1, 2] - 2 * icc[1, 0])
            except:
                dx, dy = 0, 0 #  peak is near the edge and one of the above values does not exist

            xySubShift = xySubShift - globalShift
            xyShift = xyShift + (xySubShift + np.array([dx, dy])) / upsample_factor

    return xyShift


def upsampleFFT(image_init, upsample_factor):
    """This does a Fourier upsample of the imageInit. imageInit is the Fourier transform of the correlation image.
    The function returns the real space correlation image that has been Fourier upsampled by the upsample_factor.
    An upsample factor of 2 is generally sufficient.

    The way it works is that it embeds imageInit in a larger array of zeros, then does the inverse Fourier transform to return the Fourier upsampled image in real space.

    Parameters
    ----------
        image_init : ndarray  complex
            The image to be Fourier upsampled. This should be in the Fourier domain.
        upsample_factor : int
            THe upsample factor (usually 2).

    Returns
    -------
        imageUpsampleReal : ndarray complex
            The inverse Fourier transform of imageInit upsampled by the upsampleFactor.
    OLD
    ---
    imageSize = imageInit.shape
    imageUpsample = np.zeros(tuple((i*upsampleFactor for i in imageSize))) + 0j
    imageUpsample[:imageSize[0], :imageSize[1]] = imageInit
    imageUpsample = np.roll(np.roll(imageUpsample, -int(imageSize[0]/2), 0), -int(imageSize[1]/2),1)
    imageUpsampleReal = np.real(np.fft.ifft2(imageUpsample))
    return imageUpsampleReal
    """
    
    ss = [int(ii * upsample_factor / 4) for ii in image_init.shape]  # pad size
    imageInit2 = np.pad(np.fft.fftshift(image_init), ss, mode='constant')  # pad the FFT
    image_upsample_real = np.real(np.fft.ifftn(np.fft.ifftshift(imageInit2)))  # inverse FFT

    return image_upsample_real


def dftUpsample(image_corr, upsample_factor, xy_shift):
    """
    This performs a matrix multiply DFT around a small neighboring region of the initial correlation peak.
    By using the matrix multiply DFT to do the Fourier upsampling, the efficiency is greatly improved.
    This is adapted from the subfunction dftups found in the dftregistration function on the Matlab File Exchange.

    https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

    The matrix multiplication DFT is from
    Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms,"
    Opt. Lett. 33, 156-158 (2008). http://www.sciencedirect.com/science/article/pii/S0045790612000778

    Parameters
    ----------
        image_corr : ndarray
            Correlation image between two images in Fourier space.
        upsample_factor : int
            Scalar integer of how much to upsample.
        xy_shift : list of 2 floats
            Single pixel shift between images previously computed. Used to center the matrix multiplication
            on the correlation peak.

    Returns
    -------
        image_upsample : ndarray
            Upsampled image from region around correlation peak.

    """
    imageSize = image_corr.shape
    pixelRadius = 1.5
    numRow = np.ceil(pixelRadius * upsample_factor)
    numCol = numRow

    colKern = np.exp((-1j * 2 * np.pi / (imageSize[1] * upsample_factor)) *
                     (np.fft.ifftshift( (np.arange(imageSize[1]))) -
                     np.floor(imageSize[1]/2)) *
                     (np.arange(numCol) - xy_shift[1])[:, np.newaxis]
                     )

    rowKern = np.exp(
    (-1j * 2 * np.pi / (imageSize[0] * upsample_factor))
    * (np.arange(numRow) - xy_shift[0])
    * (np.fft.ifftshift(np.arange(imageSize[0]))
    - np.floor(imageSize[0]/2))[:, np.newaxis]
    ) # Comment from above applies.

    image_upsample = np.real(np.dot(np.dot(rowKern.transpose(), image_corr), colKern.transpose()))

    return image_upsample


def imageShifter(g1, xy_shift):
    """
    Multiply im by a plane wave that has the real space effect of shifting ifft2(G2) by [x, y] pixels.

    Parameters
   -----------
        g1 : complex ndarray
            The Fourier transform of an image.
        xy_shift : list
            A two element list of the shifts along each axis.

    Returns
    -------
        G2shift : complex ndarray
            Fourier shifted image FFT

    Example
    -------
        >> shiftIm0 = np.real(np.fft.ifft2(multicorr.imageShifter(np.fft.fft2(im0),[11.1,22.2])))
        >> plt.imshow(shiftIm0)

    """
    # Check that the inputs are complex FFTs (common error)
    if not np.iscomplexobj(g1):
        raise TypeError('g1 must be complex FFTs.')

    imageSize = g1.shape
    qx = np.fft.fftfreq(imageSize[0], 1)  # does this need to be a column vector
    if imageSize[1] == imageSize[0]:
        qy = qx
    else:
        qy = np.fft.fftfreq(imageSize[1], 1)

    G2shift = np.multiply(g1, np.outer(np.exp(-2j * np.pi * qx * xy_shift[0]), np.exp(-2j * np.pi * qy * xy_shift[1])))

    return G2shift
