"""
Module to calculate radial profiles.

The **settings** dict used for the evaluation has the following entries:

    * lmax_r (int):    Radius used for local maximum detection in [px].
    * lmax_thresh (int):    Threshold for local maximum detection in [intensity].
    * lmax_cinit (int, int):    Initial guess for center in [px].
    * lmax_range (float, float):    r range to cinit used to filter local maxima [dims].
    * ns (int, ...):    Distortion orders to correct.
    * fit_rrange (float, float):    r range used to fit radial profile [dims].
    * back_xs (float, ...):    Fixpoints for fitting background [dims].
    * back_xswidth (float):    Range around fixpoints to take into account [dims].
    * back_init (float):    Initial guess for subtracting background.
    * fit_funcs (str, ...):    List of functions to model radial profile.
    * fit_init (float, ...):    Initial guess for fitting.

    optional (set to None to use defaults):

        * plt_imgminmax (float, float):    Relative range to plot the img.
        * rad_rmax (float):    Max r to be used in radial profile [dims].
        * rad_dr (float): Stepsize for r-axis in radial profile [dims].
        * rad_sigma (float): Sigma for Gaussian used as kernel density estimator [dims].
        * mask (np.ndarray): Binary image as img, 0 for pixels to exclude.
        * fit_maxfev (int): Maxfev forwarded to scipy optimize.

"""

import copy
import numpy as np
import scipy.ndimage
import scipy.interpolate

import ncempy.algo.math
import ncempy.algo.distortion
import ncempy.viz


def calc_polarcoords(center, dims, ns=None, dists=None):
    """Calculate the polar coordinates for an image of given shape.

    Center is assumed to be in real coordinates (if not just fake the dims).
    Distortions are corrected if ns and corresponding dists are given.

    Parameters
    ----------
        center : np.ndarray/tuple
            Center of polar coordinate system.
        dims : tuple
            Tuple of dimensions.
        ns : tuple
            List of distortion orders to correct for.
        dists : np.ndarray
            Parameters for distortions.

    Returns
    --------
        : tuple
            Polar coordinates (r, theta) as two np.ndarrays with same dimensions as original image.

    """

    # check input
    try:
        # check center 
        center = np.array(center)
        center = np.reshape(center, 2)

        # check if enough dims available
        assert (len(dims) >= 2)
        assert (len(dims[0]) == 3)

        # check orders
        if ns is not None:
            assert (len(ns) >= 1)

            # check dists
            assert (dists.shape[0] == len(ns) * 2 + 1)
    except TypeError:
        raise TypeError('Something wrong with the input!')

    # create coordinates
    xx, yy = np.meshgrid(dims[0][0], dims[1][0], indexing='ij')

    # calculate polar coordinate system
    rs = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    thes = np.arctan2(yy - center[1], xx - center[0])

    # correct for distortions
    if ns is not None:
        for i in range(len(ns)):
            rs /= ncempy.algo.distortion.rad_dis(thes, dists[i * 2 + 1], dists[i * 2 + 2], ns[i])

    return rs, thes


def correct_distortion(img, dims, center, ns, dists):
    """Give corrected version of diffraction pattern with respect to distortions.

    Parameters
    -----------
        img : np.ndarray
            Image (2D).
        dims : tuple
            Dimensions tuple.
        center : np.ndarray or tuple
            Center to be used.
        ns : tuple
            List of distortion orders.
        dists : np.ndarray
            Distortion parameters.

    Returns
    --------
        : np.ndarray
           Corrected diffraction pattern.

    """

    # check input
    try:
        assert (isinstance(img, np.ndarray))

        # check center 
        center = np.array(center)
        center = np.reshape(center, 2)

        # check if enough dims available
        assert (len(dims) >= 2)
        assert (len(dims[0]) == 3)

        # check orders
        assert (len(ns) >= 1)

        # check dists
        assert (dists.shape[0] == len(ns) * 2 + 1)
    except TypeError:
        raise TypeError('Something wrong with the input!')

    rs, thes = calc_polarcoords(center, dims)

    # anti distort
    for i in range(len(ns)):
        rs *= ncempy.algo.distortion.rad_dis(thes, dists[i * 2 + 1], dists[i * 2 + 2], ns[i])

    dis_xx = rs * np.cos(thes) + center[0]
    dis_yy = rs * np.sin(thes) + center[1]

    f = scipy.interpolate.RectBivariateSpline(dims[0][0], dims[1][0], img)

    return f.ev(dis_xx, dis_yy)


def calc_radialprofile(img, rs, rMax, dr, rsigma, mask=None):
    """Calculate the radial profile using Gaussian kernel density estimator.

    It is suggested to use an rMax such that all directions are still in the image, otherwise outer areas will contribute differently and lead to different signal-to-noise. A value of dr corresponding to 1/10 px and rsigma corresponding to 1 px are good parameters.

    Parameters
    -----------
        img : np.ndarray
            Image to take radial profile of intensity from.
        rs : np.ndarray
            Array containing the radial distance, same shape as img.
        rMax : float
              Maximum radial distance used for the radial profile.
        dr : float
            Step size for r-axis of radial distance.
        rsigma : float
            Sigma for Gaussian used as kernel density estimator.
        mask : np.ndarray
            Mask image, values at 0 are excluded from radial profile.

    Returns
    --------
        : tuple
            Tuple of radial and intensity axes.

    """

    # check input
    try:
        assert (isinstance(img, np.ndarray))
        assert (isinstance(rs, np.ndarray))
        assert (np.array_equal(img.shape, rs.shape))

        rMax = float(rMax)
        dr = float(dr)
        rsigma = float(rsigma)

        if not mask is None:
            assert (isinstance(mask, np.ndarray))
            assert (np.array_equal(img.shape, mask.shape))

    except:
        raise TypeError('Something wrong with input.')

    # process mask
    if mask is None:
        mask = np.ones(img.shape)
    mask = mask.astype('float64')
    mask[np.where(mask > 0)] = 1
    mask[np.where(mask == 0)] = np.nan

    # prepare radial axis for hist  
    rBins = np.arange(0, rMax, dr)

    radialIndex = np.round(rs / dr) + 1
    rBinMax = len(rBins)
    sel = (radialIndex <= rBinMax)

    sel = np.logical_and(sel, np.logical_not(np.isnan(mask)))

    signal = np.histogram(rs[sel], rBins, weights=img[sel])
    count = np.histogram(rs[sel], rBins, weights=np.ones(img[sel].shape))

    signal_sm = scipy.ndimage.gaussian_filter1d(signal[0], rsigma / dr)
    count_sm = scipy.ndimage.gaussian_filter1d(count[0], rsigma / dr)

    # masked regions lead to 0 in count_sm, divide produces nans, just ignore the warning    
    old_err_state = np.seterr(divide='ignore', invalid='ignore')
    signal_sm = np.divide(signal_sm, count_sm)
    np.seterr(**old_err_state)

    return rBins[:-1], signal_sm


def residuals_fit(param, r, intens, funcs):
    """Residual function to fit radial profile with flexibility.

    The arguments for the single functions are automatically cut from the fit parameters.

    Parameters
    -----------
        param : np.ndarray
            Fit parameters.
        r : np.ndarray
            r-axis.
        intens : np.ndarray
            Intensity-axis.
        funcs : tuple
            List of functions to include.

    Returns
    --------
        : np.ndarray
            Residuals.

    """

    return intens - ncempy.algo.math.sum_functions(r, funcs, param)


def fit_radialprofile(r, intens, funcs, init_guess, maxfev=None):
    """Fit the radial profile.

    Convenience wrapper for fitting.

    Parameters
    -----------
        r : np.ndarray
            r-axis of radial profile.
        intens : np.ndarray
            Intensity-axis of radial profile.
        funcs : tuple
            List of functions.
        init_guess : np.ndarray
            Initial guess for parameters of functions in funcs.
        maxfev : int
           maximum function calls (for scipy.optimize.leastsq)

    Returns
    --------
        : np.ndarray
            Optimized parameters.

    """

    try:
        # check data
        assert (isinstance(r, np.ndarray))
        assert (isinstance(intens, np.ndarray))
        assert (np.array_equal(r.shape, intens.shape))

        # funcs and params
        assert (len(funcs) >= 1)
        for i in range(len(funcs)):
            assert (funcs[i] in ncempy.algo.math.lkp_funcs)

        init_guess = np.array(init_guess)
        init_guess = np.reshape(init_guess, sum(map(lambda x: ncempy.algo.math.lkp_funcs[x][1], funcs)))

    except:
        raise TypeError('Something wrong with the input!')

    if maxfev is None:
        maxfev = 1000

    popt, flag = scipy.optimize.leastsq(residuals_fit, init_guess, args=(r, intens, funcs), maxfev=maxfev)

    if flag not in [1, 2, 3, 4]:
        print('WARNING: fitting of radial profile failed.')

    return popt


def run_singleImage(img, dims, settings, show=False):
    """Evaluate a single ring diffraction pattern with given settings.

    Parameters
    -----------
        img : np.ndarray
            Image.
        dims : tuple
            Corresponding dim vectors.
        settings : dict
            Dict of settings necessary for the evaluation.
        show : bool
            Set to directly show plots interactively.

    Returns
    --------
        (np.ndarray):    Optimized parameters of fitting the radial profile according to the settings.

    """

    try:
        assert (isinstance(img, np.ndarray))

        # check if dims availabel
        assert (len(dims) >= 1)
        assert (len(dims[0]) == 3)

        assert (type(settings) is dict)

    except:
        raise RuntimeError('Something wrong with the input')

    # create a copy of settings to decouple
    mysettings = copy.deepcopy(settings)

    # get local maxima an turn them into real space coords
    points = ncempy.algo.local_max.local_max(img, mysettings['lmax_r'], mysettings['lmax_thresh'])
    points = ncempy.algo.local_max.points_todim(points, dims)

    # convert center to real space
    center_init = ncempy.algo.local_max.points_todim(mysettings['lmax_cinit'], dims)

    # filter to single ring
    points = ncempy.algo.distortion.filter_ring(points, center_init, mysettings['lmax_range'])

    if show:
        if settings['plt_imgminmax'] is None:
            mysettings['plt_imgminmax'] = (0., 1.)

        plot = ncempy.viz.plot_points(img, points, vminmax=mysettings['plt_imgminmax'], dims=dims, invert=True,
                                      show=show)

    # optimize center
    center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=show)

    # fit distortions
    points_plr = ncempy.algo.distortion.points_topolar(points, center)
    dists = ncempy.algo.distortion.optimize_distortion(points_plr, mysettings['ns'])
    if show:
        plot = ncempy.viz.plot_distpolar(points_plr, dims, dists, mysettings['ns'], show=show)

    # calc coordinates in optimized system
    rs, thes = ncempy.algo.radial_profile.calc_polarcoords(center, dims, mysettings['ns'], dists)

    if settings['rad_rmax'] is None:
        mysettings['rad_rmax'] = np.abs(dims[0][0][0] - dims[0][0][1]) * np.min(img.shape) / 2.0
    if settings['rad_dr'] is None:
        mysettings['rad_dr'] = np.abs(dims[0][0][0] - dims[0][0][1]) / 10.
    if settings['rad_sigma'] is None:
        mysettings['rad_sigma'] = np.abs(dims[0][0][0] - dims[0][0][1])

    # extract radial profile
    R, I = ncempy.algo.radial_profile.calc_radialprofile(img, rs, mysettings['rad_rmax'], mysettings['rad_dr'],
                                                         mysettings['rad_sigma'], mask=mysettings['mask'])

    # cut radial profile
    sel = (R >= mysettings['fit_rrange'][0]) * (R <= mysettings['fit_rrange'][1])
    I = I[sel]
    R = R[sel]

    rawRI = np.copy(np.array([R, I]).transpose())

    # subtract a power law background fitted to specific points
    fit_R = np.array([])
    fit_I = np.array([])
    for xpoint in mysettings['back_xs']:
        ix = np.where(np.abs(R - xpoint) < mysettings['back_xswidth'])
        fit_R = np.append(fit_R, R[ix])
        fit_I = np.append(fit_I, I[ix])

    funcs_back = ['const', 'powlaw']
    res_back = ncempy.algo.radial_profile.fit_radialprofile(fit_R, fit_I, funcs_back, mysettings['back_init'],
                                                            maxfev=1000)
    if show:
        plot = ncempy.viz.plot_fit(R, I, dims, funcs_back, res_back, show=show)

    I = I - ncempy.algo.math.sum_functions(R, funcs_back, res_back)

    # fit
    res = ncempy.algo.radial_profile.fit_radialprofile(R, I, mysettings['fit_funcs'], mysettings['fit_init'],
                                                       maxfev=mysettings['fit_maxfev'])

    if show:
        plot = ncempy.viz.plot_fit(R, I, dims, mysettings['fit_funcs'], res, show=show)

    return np.array([R, I]).transpose(), res, center, dists, rawRI, res_back, mysettings
