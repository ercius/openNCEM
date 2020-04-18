import numpy as np
from scipy import ndimage


def imageCrossCor(image, reference, real_filter=1, k_filter=1):
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

    if np.iscomplex(image) or np.iscomplex(reference):
        raise TypeError("Images mst be real")

    # Calculate FFT of images
    image_f = np.fft.fft2((image - np.mean(image)) * real_filter)
    reference_f = np.fft.fft2((reference - np.mean(reference)) * real_filter)

    # Calculate cross-correlation
    xcor = abs(np.fft.fftshift(np.fft.ifft2(np.conj(image_f) * reference_f * k_filter)))

    return xcor


def imageCrossCorRealRoll(image, reference, real_filter=1, k_filter=1):
    """ Align image to reference by cross-correlation.
    Uses the real FFT for ~2x speed improvement. The k_filter must have
    a shape that matches the np.fft.rfft2() of image and reference.
    Shift the input image using np.roll to be reversible.

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

        real_filter : ndarray, optional, default = 1 (no filter)
            A real space filter applied to image and reference before
            calculating the shift.

        k_filter : ndarray, optional, default = 1 (no filter)
            A Fourier space filter applied to the fourier transform of
            image and reference before calculating the cross-correlation.

    Returns
    ------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.

    """

    image_f = np.fft.rfft2((image - np.mean(image)) * real_filter)
    reference_f = np.fft.rfft2((reference - np.mean(reference)) * real_filter)
    xcor = abs(np.fft.irfft2(np.conj(image_f) * reference_f * k_filter))

    shifts = np.unravel_index(xcor.argmax(), xcor.shape)
    shifts = [int(i) for i in shifts]  # convert to integers

    # shift image using roll to be reversible
    output = np.roll(image, shifts[0], axis=0)
    output = np.roll(output, shifts[1], axis=1)

    return output, shifts


def imageCrossCorRealShift(image, reference, real_filter=1, k_filter=1, verbose=False):
    """ Align image to reference by cross-correlation.
    Uses the real FFT for ~2x speed improvement. The k_filter must have
    a shape that matches the np.fft.rfft2() of image and reference.
    Uses scipy.ndimage.shift() to shift the image and remove border pixels.

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

        verbose : bool
            Plots the cross-correlation using matplotlib for debugging purposes.

    Returns
    ------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.
    """
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

    # shift image using ndimage.shift
    output = ndimage.interpolation.shift(image, shifts, order=0)

    return output, shifts


def imageCrossPhaseRealShift(image, reference, real_filter=1, k_filter=1, verbose=False):
    """ Align image to reference by phase-correlation.
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

        verbose : bool
            Plots the cross-correlation using matplotlib for debugging purposes.

    Returns
    ------
        : tuple, (ndarray, tuple)
            A tuple containing the shifted image and the shifts applied.


    """

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

    # shift image using ndimage's shift
    output = ndimage.interpolation.shift(image, shifts, order=0)

    return output, shifts


def stackCrossCorReal(tilt_series, angles, real_filter=1, k_filter=1, use_center_ref=True):
    """ Align a tilt series of images by cross-correlation. The
    function uses the min(abs(angles) to determine the angle closest
    to zero degree tilt. All images are aliged to this start image from
    start -> max(angles) and then start -> min(angles). The function
    uses np.fft.rfft2() to increase the speed.

    Adapted from tomviz v1.0.0.

    Parameters
    ----------
        tilt_series : ndarray, 3D
            The stack of images to align. Shape [num, Y, X]

        angles : ndarray, 1D
            The angles in degrees for each image in the stack.
            Shape should be [num,].

        real_filter : ndarray, optional, default = 1
            A real space filter to apply before cross-correlation
            of each image. Shape must be [Y, X]

        k_filter : ndarray, optional, default = 1
            A Fourier space filter to apply before cross-correlation.
            Shape must be [Y, X/2 + 1]

        use_center_ref : bool
            Use the central image as the reference for all images

    Returns
    -------
        : tuple
            A tuple containing the aligned images as a 3D ndarray of shape
            [num, Y, X] and shifts as a 2D ndarray of shape [num, 2]
    """

    # Determine reference image index
    zeroDegreeTiltImage = np.where(angles == 0)
    if zeroDegreeTiltImage[0]:
        referenceIndex = zeroDegreeTiltImage[0]
    else:
        print('No zero degree image. Using middle image.')
        referenceIndex = tilt_series.shape[0] // 2

    # Pre-allocate the arrays
    aligned = np.zeros_like(tilt_series)  # shifted data array
    shifts = np.zeros((tilt_series.shape[0], 2))  # the applied shifts

    # save reference image
    aligned[referenceIndex, :, :] = tilt_series[referenceIndex, :, :]

    # Align positive angles
    j = referenceIndex
    for i in range(referenceIndex, tilt_series.shape[0] - 1):
        if not use_center_ref:
            j = i
        output, sh = imageCrossCorRealShift(tilt_series[i + 1, :, :], tilt_series[j, :, :], real_filter, k_filter)
        aligned[i + 1, :, :] = output
        shifts[i + 1, :] = sh

    # Align negative angles
    for i in range(referenceIndex, 0, -1):
        if not use_center_ref:
            j = i
        output, sh = imageCrossCorRealShift(tilt_series[i - 1, :, :], tilt_series[j, :, :], real_filter, k_filter)
        aligned[i - 1, :, :] = output
        shifts[i - 1, :] = sh

    return aligned, shifts
