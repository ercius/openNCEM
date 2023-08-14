import numpy as np

from ncempy.eval import multicorr


def stack_align(stack, align_type='static', **kwargs):
    """ Align a series of images by correlation using multicorr. All images are aligned to the start image or as
    dynamically starting with the start image and then each successive image. keyword arguments are passed to
    multicorr. Shifting is done in Fourier space which is very accurate, but wraps edges.

    Parameters
    ----------
        stack : ndarray, 3D,
            Stack of images to align. Shape [num, Y, X]

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

    # Align positive angles
    ref_fft = np.fft.fft2(stack[0, :, :])
    ref_sh = np.zeros((2,))
    for ii in range(1, stack.shape[0]):
        cur_fft = np.fft.fft2(stack[ii, :, :])
        sh = multicorr.multicorr(cur_fft, ref_fft, **kwargs)
        sh += ref_sh
        image_shifted = np.real(np.fft.ifft2(multicorr.imageShifter(cur_fft, sh)))
        aligned[ii, :, :] = image_shifted
        shifts[ii, :] = sh

        if align_type == 'dynamic':
            ref_fft = cur_fft
            ref_sh = sh

    return aligned, shifts
