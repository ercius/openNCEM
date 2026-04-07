# Convenience EELS background fitting and analysis functions. Hyperspy and other dedicated packages are much more
# powerful, but simple functions like these here can be used in cases where more simple or more complex analysis
# can benefit from access to low level functionality. For example, a Gaussian-Lorentz background fit for low-loss
# EELS is included which allows for pre- and post-edge values to be used in the fit.

import numpy as np
from scipy import optimize as opt

from .gaussND import gauss1D
from .gaussND import lorentz1D


def powerFit(sig, box):
    """Fit a power low function to a signal inside a designated range.

    Parameters
    ----------
    sig: ndarray, 1D
        The signal to fit the power law to. Usually the spectral intensity.
    box: 2-tuple
        The range to fit of the form (low, high) where low and high are in terms of the lowest and highest pixel number
        to include in the range to fit.

    Returns
    -------
    : tuple
        A tuple is returned. The first element contains the values of the power law background fit for all pixel values above the lowest pixel
        in the range and the end of the spectrum. The second element are the pixel positions of the fit.
    """
    eLoss = np.linspace(box[0], box[1] - 1, int(np.abs(box[1] - box[0])))
    try:
        yy = sig.copy()
        yy[yy < 0] = 0
        yy += 1

        # TODO Replace polyfit with numpy.Polynomial.fit
        pp = np.polyfit(np.log(eLoss), np.log(yy[box[0]:box[1]]), 1)  # log-log fit
    except:
        print('fitting problem')
        raise

    pixel_range = np.linspace(box[0], sig.shape[0], np.int(np.abs(sig.shape[0] - box[0])))
    try:
        background = np.exp(pp[1]) * pixel_range ** pp[0]
    except:
        print('background creation problem')
        background = np.zeros_like(pixel_range)
    return background, pixel_range


def pre_post_fit(energy_axis, spectra, pre, post, initial_parameters=(0, .15, .15, 0.5, 0)):
    """ Fit a Gaussian-Lorentz to pre- and post-edge signals. This is useful for low loss EELS

    The function is of the form:
    GL = n * Gaussian(x0, sigma) + (1-n) * Lorentz(x0, w) + C

    Parameters
    ----------
    energy_axis: ndarray, 1D
        The energy loss values for each spectral intensity.
    spectra: ndarray, 1D
        The values of the spectra at the corresponding energy loss values.
    pre: 2-tuple
        Start and end pixels for the pre-edge fitting range (start, end)
    post: 2-tuple
        Start and end pixels for the post-edge fitting range (start, end)
    initial_parameters: tuple
        The initial guess for the parameters of the form (x0, sigma, w, n, C) as described in the function definition.
        All values are in eV except n which is dimensionless but should be [0, 1]. Default initial parameters are for
        a zero loss peak centered at 0 eV with a width of 0.15 eV, 50% mixing of the gaussian/lorentz and zero offset.

    Returns
    -------
    : 2-tuple
        A tuple contains two 1D ndarrays. The first is the energy loss values for the signal and the second is the
        background subtracted signal.

    """
    # Put the pre and post signals in the same arrays
    eLoss_fit = np.concatenate((energy_axis[pre[0]:pre[1]], energy_axis[post[0]:post[1]]))
    sig_fit = np.concatenate((spectra[pre[0]:pre[1]], spectra[post[0]:post[1]]))

    # Gaussian-Lorentz fit to the pre- and post-edge
    # TODO allow for other functions
    GL = lambda x, x0, sigma, w, n, C: n * gauss1D(x, x0, sigma) + (1 - n) * lorentz1D(x, x0, w) + C
    pp, r = opt.curve_fit(GL, eLoss_fit, sig_fit, p0=initial_parameters, maxfev=20000)

    # Create background using fitted values
    bgnd_eloss = np.linspace(energy_axis[pre[0]], energy_axis[post[1]], energy_axis[pre[0]:post[1]].shape[0])
    bgnd = GL(bgnd_eloss, *pp)

    # Subtract the background from the signal
    sig_subtracted = spectra[pre[0]:post[1]] - bgnd
    return bgnd_eloss, sig_subtracted


def adaptive_background():
    """An adaptive background which sum pixels in a neighborhood of each probe position to improve SNR for the
    background fit. This should only be used when the sampling is higher than the variation in the background signal."""
    pass


def map_background(spectra, bgnd_box):
    """ Power law background subtraction to a set of spectra.

    Parameters
    ----------
    spectra : ndarray, 2D or 3D
        Spectra as a 2D array (line scan) with shape (N, energy_loss) or a 3D array (spectrum image) with shape
        (N, M, energy_loss).

    bgnd_box : 2-tuple of ints
        The starting and ending point (in pixels) for the powerLaw background subtraction.

    Returns
    -------
    : ndarray, 2D or 3D
        The background subtracted spectra.
    """
    if spectra.ndim == 3:
        mm = spectra.reshape((spectra.shape[0] * spectra.shape[1], spectra.shape[2]))
    elif spectra.ndim == 2:
        mm = spectra
    else:
        print('Incorrect spectra size')
        return
    mm_sub = np.zeros_like(mm)
    for ii, sig in enumerate(mm):
        bgnd, eLoss0 = powerFit(sig, bgnd_box)
        mm_sub[ii, bgnd_box[0]:] = sig[bgnd_box[0]:] - bgnd

    if spectra.ndim == 3:
        mm_sub = mm_sub.reshape(*spectra.shape)

    return mm_sub


def map_sum(spectra, sig_box):
    """ Sum the signal in the spectrum between two pixel values. This will essentially produce a compositional map of
    the element of interset.

    Parameters
    ----------
    spectra : ndarray, 2D or 3D
        Spectra as a 2D array (line scan) with shape (N, energy_loss) or a 3D array (spectrum image) with shape
        (N, M, energy_loss).
    sig_box : 2-tuple of ints
        The starting and ending point (in pixels) for the signal summation.

    Returns
    -------
    : ndarray, 1D or 2D
        For a line scan the return is a 1D set of points and for a spectrum image the returned 2D array is a map.


    """
    sig_map = np.sum(spectra[..., sig_box[0]:sig_box[1]], axis=-1)
    return sig_map


def create_rgb(maps, normalize=True):
    """ Create an image that combines the input maps into a color image. Each map will be used as the red, green, and
    blue channels of the final image.

    Parameters
    ----------
    maps : ndarray, 3D
        A 3D array of shape (3, M, N) where the first two dimensions are the color maps and the thrid dimension
        indicates the red, green, and blue channels in that order.
    normalize : bool
        Set the maximum pixel of each image to 1 for each channel independently. Default is True.

    Returns
    -------
    : ndarray, 3D
        A 3D array of float values. If you pass this to matplotlib.pyplot.imshow or write it to a PNG file it will be
        displayed as a color image where the red, green, blue channels are properly mixed.
    """
    if maps.ndim != 3:
        raise ValueError("map array must be 3-dimensional.")

    rgb = np.zeros((maps[0].shape[0], maps[0].shape[1], 3), dtype=np.float32)
    for ii, m in enumerate(maps):
        if normalize:
            m0 = m - m.min()
            m0 = m0 / m0.max()
        else:
            m0 = m
        rgb[:, :, ii] = m0
    return rgb
