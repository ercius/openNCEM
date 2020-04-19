import numpy as np

import ncempy.algo.multicorr as multicorr


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

    if align_type is not 'static' and align_type is not 'dynamic':
        raise KeyError('Incorrect align type. Must be static or dynamic')

    # Pre-allocate the arrays
    aligned = np.zeros_like(stack)  # shifted data array
    shifts = np.zeros((stack.shape[0], 2))  # the applied shifts

    aligned[0, :, :] = stack[0, :, :]

    # Align positive angles
    ref_fft = np.fft.fft2(stack[0, :, :])
    for ii in range(1, stack.shape[0]):
        cur_fft = np.fft.fft2(stack[ii, :, :])
        sh = multicorr.multicorr(cur_fft, ref_fft, **kwargs)
        image_shifted = np.real(np.fft.ifft2(multicorr.imageShifter(cur_fft, sh)))
        aligned[ii, :, :] = image_shifted
        shifts[ii, :] = sh

        if align_type is 'dynamic':
            jj = ii - 1
            ref_fft = cur_fft

    return aligned, shifts


if __name__ == '__main__':

    import ncempy.io as nio
    from scipy import ndimage

    up = 2
    sh0 = (4.5, 5.5)
    method0 = 'hybrid'

    print('Upsample factor = {}'.format(up))
    print('Applied shift = {}'.format(sh0))

    with nio.emd.fileEMD('C:/Users/linol/Data/Acquisition_18.emd') as f0:
        dd, md = f0.get_emdgroup(f0.list_emds[0])

    stack0 = np.zeros((5, *dd.shape), dtype=dd.dtype)
    shifts0 = np.asarray(((range(4, -6, -2)), (range(6, -8, -3)))).T
    for ii, s in enumerate(shifts0):
        stack0[ii, :, :] = ndimage.shift(dd, s, mode='mirror')
    out_stack, out_stack_shifts = stack_align(stack0, align_type='static',
                                              upsample_factor=up, method='hybrid')
    print('shifts0 = {}'.format(shifts0))
    print('stack shifts = {}'.format(out_stack_shifts))