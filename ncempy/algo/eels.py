# EELS background fitting and analysis functions
import numpy as np


def powerFit(sig, box):
    """A powerlaw fit"""
    eLoss = np.linspace(box[0], box[1] - 1, int(np.abs(box[1] - box[0])))
    try:
        yy = sig  # .copy()#[int(box.left()):int(box.right())]
        yy[yy < 0] = 0
        yy += 1

        pp = np.polyfit(np.log(eLoss), np.log(yy[box[0]:box[1]]), 1)  # log-log fit
    except:
        print('fitting problem')
        pp = (0, 0)
        raise
    fullEloss = np.linspace(box[0], sig.shape[0], np.int(np.abs(sig.shape[0] - box[0])))
    try:
        bgnd = np.exp(pp[1]) * fullEloss ** pp[0]
    except:
        print('bgnd creation problem')
        bgnd = np.zeros_like(fullEloss)
    return bgnd, fullEloss


def gauss1D(x, x0, sigma):
    """ Returns the value of a gaussian at a 2D set of points for the given standard deviation with maximum
    normalized to 1.

    Parameters
    ----------
        x: ndarray
            A vector of size (N,) of points for the x values
        x0: float
            The center of the Gaussian.
        sigma: float
            Standard deviation of the Gaussian.

    Returns
    -------
        g: ndarray
            A vector of size (N,) of the Guassian distribution evaluated
            at the input x values.

    Note
    ----
        Calculate the Half width at half maximum as
        HWHM = sqrt(2*log(2))*stDev ~ 1.18*stDev or
        0.5*size(x)*stDev if x goes from -1 to 1

    """
    return np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def lorentz1D(x, x0, w):
    """ Returns the probability density function of the Lorentzian (aka Cauchy) function with maximum at 1.

    Parameters
    ----------
        x: ndarray or list
            A 1D vector of points for the dependent variable
        x0: float
            The center of the peak maximum
        w: float
            The parameter that modifies the width of the distribution.

    Returns
    -------
        l: ndarray
            A vector of size (N,) of the Lorentzian distribution evaluated
            at the input x values.

    """
    return w ** 2 / ((x - x0) ** 2 + (w) ** 2)


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