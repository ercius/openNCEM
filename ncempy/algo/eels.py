# EELS background fitting and analysis functions
import numpy as np
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

def pre_post_fit(energy_axis, spectra, pre, post, initial_parameters):
    """Fit to pre and post edge. Most useufl for low loss EELS"""
    # Put the pre and post signals in the same arrays
    eLoss_fit = np.concatenate((energy_axis[pre[0]:pre[1]], energy_axis[post[0]:post[1]]))
    sig_fit = np.concatenate((spectra[pre[0]:pre[1]], spectra[post[0]:post[1]]))

    # Powerlaw fit to the pre and postedge
    GL = lambda x, x0, sigma, w, n, C: n * gauss1D(x, x0, sigma) + (1 - n) * lorentz1D(x, x0, w) + C
    pp, r = opt.curve_fit(GL, eLoss_fit, sig_fit, p0=initial_parameters, maxfev=20000)

    # Create background using fit
    bgnd_eloss = np.linspace(energy_axis[pre[0]], energy_axis[post[1]], energy_axis[pre[0]:post[1]].shape[0])
    bgnd = GL(bgnd_eloss, *pp)

    # Subtract the background from the signal
    sig_subtracted = spectra[pre[0]:post[1]] - bgnd
    return bgnd_eloss, sig_subtracted