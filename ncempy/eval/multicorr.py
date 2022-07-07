import numpy as np
from ..algo.multicorr_funcs import *


def multicorr(g1, g2, method='cross', upsample_factor=1, verbose=False):
    """Align a reference to an image by cross correlation. The template
    and the image must have the same size.

    The function takes in FFTs so that any FFT algorithm can be used to
    transform the image and template (fft2, mkl, scipack, etc.)

    Parameters
    ----------
        g1 : complex ndarray
            Fourier transform of reference image.
        g2 : complex ndarray
            Fourier transform of the image to register (the kernel).
        method : str, optional
            The correlation method to use. Must be 'phase' or 'cross' or 'hybrid' (default = 'cross')
        upsample_factor : int
            Upsample factor for subpixel precision of cross correlation. (default = 1)
        verbose : bool, default is False
            Print output.

    Returns
    -------
        xyShift : list of floats
            The shift between G1 and G2 in pixels.

    Example
    -------
        Cross correlate two images already stored as ndarrays. You must input the FFT
        of the images.

        >> import ncempy.algo as neval
        >> import numpy as np
        >> im0FFT = np.fft.fft2(im0)
        >> im1FFT = np.fft.fft2(im1)
        >> shifts = neval.multicorr(im0FFT, im1FFT)
    """

    # Check to make sure both G1 and G2 are arrays
    if type(g1) is not np.ndarray:
        raise TypeError('G1 must be an ndarray')
    elif type(g2) is not np.ndarray:
        raise TypeError('G2 must be an ndarray')

    # Check that the inputs are complex FFTs (common error)
    if not np.iscomplexobj(g1) or not np.iscomplexobj(g2):
        raise TypeError('G1 and G2 must be complex FFTs.')

    # Check to make sure method and upsample factor are the correct values
    if method not in ['phase', 'cross', 'hybrid']:
        print('Unknown method used, setting to cross.')
        method = 'cross'

    if type(upsample_factor) is not int and type(upsample_factor) is not float:
        print('Upsample factor is not an integer or float, setting to 1')
        upsample_factor = 1
    elif type(upsample_factor) is not int:
        print('Upsample factor is not an integer, rounding down')
        upsample_factor = int(upsample_factor)
        if upsample_factor < 1:
            print('Upsample factor is < 1, setting to 1')
            upsample_factor = 1

    if upsample_factor < 1:
        raise ValueError('upsample_factor must be >= 1')

    if verbose:
        print('upsample factor = {}'.format(upsample_factor))

    # Verify images are the same size.
    if g1.shape != g2.shape:
        raise TypeError('G1 and G2 are not the same size, G1 is {0} and G2 is {1}'.format(g1.shape, g2.shape))

    imageCorr = initial_correlation_image(g1, g2, method, verbose=verbose)
    xyShift = upsampled_correlation(imageCorr, upsample_factor, verbose=verbose)

    return xyShift
