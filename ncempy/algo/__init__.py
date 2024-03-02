
import numpy as np
from scipy import ndimage

from . import radial_profile
from . import peak_find

def rebin(im, f, funcType='sum'):
    """ Rebin a 2D array. From stackoverflow:
    https://stackoverflow.com/questions/4624112/grouping-2d-numpy-array-in-average

    This should be used carefully. There is not antialiasing applied which could
    produce odd effects for some data sets (such as high resolution data)

    Parameters
    ----------
    im : ndarray
        2D array to reduce
    f : int
        The factor to rebin by. Must be integer
    funcType : str
        The type of reduction. mean or sum are implemented.

    Returns
    -------
    : ndarray
        The 2D array with the new size.
    """

    nbig = im.shape
    nsmall = [ii // f for ii in im.shape]

    # Reshape the array so that the required neighborhood are in arrays along certain axes
    # Then average or sum those neighborhoods
    im_reshape = im.reshape([nsmall[0], nbig[0] // nsmall[0], nsmall[1], nbig[1] // nsmall[1]])
    # Reduce using different types of functions
    if funcType == 'mean':
        small = im_reshape.mean(3).mean(1)
    elif funcType == 'sum':
        small = im_reshape.sum(3).sum(1)
    # elif funcType == 'median':
    #    small = im_reshape.median(3).median(1)
    else:
        print('defaulting to sum')
        small = im_reshape.sum(3).sum(1)

    return small

##
# Image alignment codes
##

def image_cross_corr(image, reference, real_filter=1, k_filter=1):
    """ Calculate image cross-correlation. See imageCrossCorRealShift and other
    similar functions to calculate the shift and apply the shift.

    Parameters
    ----------
        image : ndarray
            The source image to align to. Should be even sized.

        reference : ndarray
            The reference image to align to the source image. Should be even sized.

        real_filter : ndarray, optional, default 1
            Filter to apply to each image in real space before FFT.

        k_filter : ndarray, optional default 1
            Filter to apply to each image in FFT space

    Returns
    -------
        : ndarray
            Cross correlation of image and reference.

    """

    if type(image) is not np.ndarray or type(reference) is not np.ndarray:
        raise TypeError("Must use ndarrays")

    if np.iscomplexobj(image) or np.iscomplexobj(reference):
        raise TypeError("Images mst be real")

    # Calculate FFT of images
    image_f = np.fft.fft2((image - np.mean(image)) * real_filter)
    reference_f = np.fft.fft2((reference - np.mean(reference)) * real_filter)

    # Calculate cross-correlation
    xcor = abs(np.fft.fftshift(np.fft.ifft2(np.conj(image_f) * reference_f * k_filter)))

    return xcor


def image_correlate(image, reference, real_filter=1, k_filter=1, shift_func='shift', verbose=False):
    """ Align image to reference by cross-correlation. Outputs shifts and shifted images.
    Uses the real FFT for ~2x speed improvement. The k_filter must have
    a shape that matches the np.fft.rfft2() of image and reference.
    Uses scipy.ndimage.shift() or np.roll to move the image. Use 'roll' to avoid losing
    data off the edge for multiple shifting operations. Use shift to avoid wrap around problems and when there
    is only one shifting operation.

    Note
    ----
        image, reference and real_filter must all have the same shape (N, M).
        k_filter must have a shape that matches the np.fft.rfft2() of
        the other inputs: (N, M/2+1)

    Parameters
    ----------
        image : ndarray
            A image as a 2D ndarray.

        reference : ndarray
            The reference image to align to.

        real_filter : ndarray, optional, default = 1
            A real space filter applied to image and reference before
            calculating the shift.

        k_filter : ndarray, optional, default = 1
            A Fourier space filter applied to the fourier transform of
            image and reference before calculating the cross-correlation.

        shift_func : str, default is 'shift'
            The function to use to shift the images. 'roll' uses np.roll and 'shift' uses ndimage.shift.

        verbose : bool
            Plots the cross-correlation using matplotlib for debugging purposes.

    Returns
    -------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.
    """
    output = None

    if shift_func != 'shift' and shift_func != 'roll':
        raise KeyError('Shift function has to be either shift or roll')

    image_f = np.fft.rfft2((image - np.mean(image)) * real_filter)
    reference_f = np.fft.rfft2((reference - np.mean(reference)) * real_filter)

    xcor = abs(np.fft.irfft2(np.conj(image_f) * reference_f * k_filter))
    if verbose:
        import matplotlib.pyplot as plt
        plt.imshow(np.fft.fftshift(xcor))
        plt.title('imageCrossCorRealShift xcor')
    shifts = np.unravel_index(np.fft.fftshift(xcor).argmax(), xcor.shape)
    shifts = (shifts[0] - xcor.shape[0] / 2, shifts[1] - xcor.shape[1] / 2)
    shifts = [int(i) for i in shifts]  # convert to integers

    if shift_func == 'shift':
        # shift image using ndimage.shift
        output = ndimage.shift(image, shifts, order=0)
    elif shift_func == 'roll':
        # shift image using roll to be reversible
        output = np.roll(image, shifts[0], axis=0)
        output = np.roll(output, shifts[1], axis=1)

    return output, shifts


def image_phase_correlate(image, reference, real_filter=1, k_filter=1, shift_func='shift', verbose=False):
    """ Align image to reference by phase-correlation. Outputs shifted images and shift.
    Uses np.fft.rfft2 for ~2x speed improvement.
    Uses scipy.ndimage.shift() to shift the image and remove border pixels.

    NOT WORKING OR TESTED YET

    Note
    ----
        image, reference and real_filter must all have the same shape (N, M).
        k_filter must have a shape that matches the np.fft.rfft2() of
        the other inputs: (N, M/2+1)

    Parameters
    ----------
        image : ndarray
            A image as a 2D ndarray.

        reference : ndarray
            The reference image to align to.

        real_filter : ndarray, optional, default = 1
            A real space filter applied to image and reference before
            calculating the shift.

        k_filter : ndarray, optional, default = 1
            A Fourier space filter applied to the fourier transform of
            image and reference before calculating the cross-correlation.

        shift_func : str, default is 'shift'
            The function to use to shift the images. 'roll' uses np.roll and 'shift' uses ndimage.shift.


        verbose : bool
            Plots the cross-correlation using matplotlib for debugging purposes.

    Returns
    ------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.

    """
    output = None

    image_f = np.fft.rfft2((image - np.mean(image)) * real_filter)
    reference_f = np.fft.rfft2((reference - np.mean(reference)) * real_filter)

    xcor = abs(np.fft.irfft2(np.conj(image_f) * reference_f * k_filter))
    pcor = xcor / (np.abs(xcor) + 0.001)
    if verbose:
        import matplotlib.pyplot as plt
        plt.imshow(np.fft.fftshift(pcor))
        plt.title('imageCrossPhaseRealShift pcor')
    shifts = np.unravel_index(np.fft.fftshift(pcor).argmax(), pcor.shape)
    shifts = (shifts[0] - pcor.shape[0] / 2, shifts[1] - pcor.shape[1] / 2)
    shifts = [int(i) for i in shifts]  # convert to integers

    if shift_func == 'shift':
        # shift image using ndimage.shift
        output = ndimage.interpolation.shift(image, shifts, order=0)
    elif shift_func == 'roll':
        # shift image using roll to be reversible
        output = np.roll(image, shifts[0], axis=0)
        output = np.roll(output, shifts[1], axis=1)

    return output, shifts


def stack_align(stack, align_type='static', real_filter=1, k_filter=1, shift_func='shift'):
    """ Align a series of images by cross-correlation. All images are aligned to the first image Uses image_correlate
    which is based on simple cross correlation.

    Notes
    -----
        You should probably use ncempy.algo.stack_align since it uses mutlicorr and is more
        functional.

    Parameters
    ----------
        stack : ndarray, 3D
            The stack of images to align. Shape [num, Y, X]

        real_filter : ndarray, optional, default = 1
            A real space filter to apply before cross-correlation
            of each image. Shape must be [Y, X]

        k_filter : ndarray, optional, default = 1
            A Fourier space filter to apply before cross-correlation.
            Shape must be [Y, X/2 + 1]

        shift_func : str, default is 'shift'
            The function to use to shift the images. 'roll' uses np.roll and 'shift' uses ndimage.shift.

        align_type: str
            static or dynamic alignment. Static aligns all images to the first image. Dynamic aligns
            each image to the previous image starting with the first image

    Returns
    -------
        : tuple, aligned stack, shifts
            A tuple containing the aligned images as a 3D ndarray of shape
            [num, Y, X] and shifts as a 2D ndarray of shape [num, 2]
    """

    if align_type != 'static' and align_type != 'dynamic':
        raise KeyError('Incorrect align type. Must be static or dynamic')

    # Pre-allocate the arrays
    aligned = np.zeros_like(stack)  # shifted data array
    shifts = np.zeros((stack.shape[0], 2))  # the applied shifts

    aligned[0, :, :] = stack[0, :, :]

    jj = 0
    ref_sh = np.zeros((2,))
    for ii in range(1, stack.shape[0]):
        output, sh = image_correlate(stack[ii, :, :], stack[jj, :, :], real_filter, k_filter, shift_func=shift_func)
        sh += ref_sh
        aligned[ii, :, :] = output
        shifts[ii, :] = sh
        if align_type == 'dynamic':
            ref_sh = sh
            jj = ii

    return aligned, shifts

def moments(im, order=3):
    """ Calculate the image raw moments.

    Parameters
    ----------
        im : ndarray
            A 2D array of the image to calculate the moments for

        order : int, default=3
            The order to calculate the moments to.

    Returns
    -------
        : ndarray
            An ndarray of (order+1, order+1)
    """
    order = int(order) + 1
    XX, YY = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    mm = np.zeros((order, order))
    for ii in range(order):
        XXii = XX ** ii
        for jj in range(order):
            mm[ii, jj] = np.sum(XXii * YY ** jj * im.astype(np.float64))
    return mm


def moments_central(im, cent=None, order=3):
    """ Calculate the image central moments.

    Parameters
    ----------
        im : ndarray
            A 2D array of the image to calculate the central moments for

        cent : tuple, 2 values, optional
            The centroid calculated from the image raw moments. If it
            is not supplied then it is calculated by moments()

        order : int, default = 3
            The order to calculate the orders to.

    Returns
    -------
        : ndarray
            An ndarray of size (order + 1, order +1)

    """
    if not cent:
        m = moments(im, order)
        cent = centroid(m)

    order = int(order) + 1

    XX, YY = np.mgrid[-cent[0]:im.shape[0] - cent[0],
                      - cent[1]:im.shape[1] - cent[1]]
    XX = XX.astype('float64')
    YY = YY.astype('float64')

    mc = np.zeros((order, order))
    for ii in range(order):
        XXii = XX ** ii
        for jj in range(order):
            mc[ii, jj] = np.sum(XXii * YY ** jj * im.astype(np.float64))
    return mc


def centroid(m):
    """ Find the centroid based on the raw moments.

    Parameters
    ----------
        m : ndarray
            The output from moments()

    Returns
    -------
        : tuple, 2 elements (X,Y)
            The position of the centroid as a 2-tuple
    """
    return m[1, 0] / m[0, 0], m[0, 1] / m[0, 0]


def moment_angle(mc):
    """ Calculate the orientation of the object from its central moments

    Parameters
    ----------
        mc : ndarray
            The array of the central moments. Output from moments_central

    Returns
    -------
        : float
            The angle in radians.

    """
    mu_20 = mc[2, 0] / mc[0, 0]  # - centroid[0]**2
    mu_02 = mc[0, 2] / mc[0, 0]  # - centroid[1]**2
    mu_11 = mc[1, 1] / mc[0, 0]  # - centroid[0]*centroid[1]

    th = 0.5 * np.arctan2(2 * mu_11, mu_20 - mu_02)
    return th

def shearImage(im, dim, shear_factor):
    """ Exact shear of an image using the frequency space of the FFT
    of the image.

    Currently only works for square images.

    Parameters
    ----------
        im: ndarray, 2D
            The image to shear.
        dim: int
            The axis to shear along
        shear_factor: float
            The amount of shear.

    Returns
    -------
        : ndarray
            The sheared image.
    """

    if im.shape[0] != im.shape[1]:
        print('Input image must be square.')
        return

    fftF = np.fft.fftshift(np.fft.fftfreq(im.shape[dim], d=1 / im.shape[dim]))
    yy, xx = np.meshgrid(fftF, fftF, indexing='ij')
    f = np.fft.fftshift(np.fft.fftn(im, axes=(dim,)), axes=(dim,))
    ff = f * np.exp(-2.0j * np.pi / im.shape[dim] * xx * yy * shear_factor)

    sh = np.fft.ifftn(np.fft.ifftshift(ff, axes=(dim,)), axes=(dim,))
    sh = np.real(sh)
    return sh


def shiftImage(im, shift):
    """ Exact shear of an image using the frequency space of the FFT
    of the image.

    Parameters
    ----------
        im: ndarray, 2D
            The image to shear.
        shift: tuple of floats
            The number of pixels to shift in each direction
    Returns
    -------
        : ndarray
            The shifted image.
    """

    fftF0 = np.fft.fftshift(np.fft.fftfreq(im.shape[0], d=1 / im.shape[0]))
    fftF1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d=1 / im.shape[1]))
    xx, yy = np.meshgrid(fftF0, fftF1, indexing='ij')
    f = np.fft.fftshift(np.fft.fftn(im))
    ff = f * np.exp(-2.0j * np.pi / im.shape[0] * xx * shift[0] +
                    -2.0j * np.pi / im.shape[1] * yy * shift[1])

    sh_im = np.fft.ifftn(np.fft.ifftshift(ff))
    sh_im = np.real(sh_im)
    return sh_im


def rotateImage(im, theta, pad=False):
    """ Use three shears in Fourier space to exactly (and reversibly) rotate
    an image.

    Currently only works for square images.

    Parameters
    ----------
        im: ndarray, 2D
            The image to rotate

        theta: float
            The angle to rotate by in radians

        pad: bool
            Add padding to the image before rotating

    """

    if im.shape[0] != im.shape[1]:
        print('Input image must be square.')
        return

    alpha = gamma = - np.tan(theta/2.)
    beta = np.sin(theta)

    if pad:
        p = im.shape[0] * np.sin(theta)
        im = np.pad(im, int(p//2+1))

    im1 = shearImage(im,  0, alpha)
    im2 = shearImage(im1, 1, beta)
    im3 = shearImage(im2, 0, gamma)
    im3 = np.real(im3)
    return im3
