import numpy as np
from scipy import ndimage


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
    ------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.
    """
    output = None

    if shift_func is not 'shift' and shift_func is not 'roll':
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
        output = ndimage.interpolation.shift(image, shifts, order=0)
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

    if align_type is not 'static' and align_type is not 'dynamic':
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
        if align_type is 'dynamic':
            ref_sh = sh
            jj = ii

    return aligned, shifts
