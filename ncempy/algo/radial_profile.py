'''
Module to calculate radial profiles.
'''

import copy
import numpy as np
import scipy.ndimage.filters
import scipy.interpolate
import matplotlib.pyplot as plt

import ncempy.algo.math
import ncempy.algo.distortion

def calc_polarcoords ( center, dims, ns=None, dists=None ):
    '''
    Calculate the polar coordinates for an image of given shape.
    
    Center is assumed to be in real coordinates (if not just fake the dims).
    Distortions are corrected if ns and corresponding dists are given.
    
    input:
    - center
    - dims
    - ns
    - dists
    
    return:
    - rs, thes
    '''
    
    # check input
    try:
        # check center 
        center = np.array(center)
        center = np.reshape(center, 2)
        
        # check if enough dims availabel
        assert(len(dims)>=2)
        assert(len(dims[0])==3)
        
        # check orders
        if not ns is None:
            assert(len(ns)>=1)
        
            # check dists
            assert(dists.shape[0] == len(ns)*2+1)
            
    except:
        raise TypeError('Something wrong with the input!')
    
    # create coordinates
    xx, yy = np.meshgrid( dims[0][0], dims[1][0], indexing='ij' )
    
    # calculate polar coordinate system
    rs = np.sqrt( np.square(xx-center[0]) + np.square(yy-center[1]) )
    thes = np.arctan2(yy-center[1], xx-center[0])
    
    # correct for distortions
    if not ns is None:
        for i in range(len(ns)):
            rs /= ncempy.algo.distortion.rad_dis(thes, dists[i*2+1], dists[i*2+2], ns[i])    
    
    return rs, thes
    
      
def correct_distortion( img, dims, center, ns, dists ):
    '''
    Give corrected version of img with respect to distortions.
    
    input:
    - img       2D np.ndarray holding image
    - dims      dimensions
    - center    center to be used
    - ns        list of distortion orders
    - dists     distortion parameters
    '''
    
    # check input
    try:
        assert(isinstance(img, np.ndarray))
    
        # check center 
        center = np.array(center)
        center = np.reshape(center, 2)
        
        # check if enough dims availabel
        assert(len(dims)>=2)
        assert(len(dims[0])==3)
        
        # check orders
        assert(len(ns)>=1)
        
        # check dists
        assert(dists.shape[0] == len(ns)*2+1)
            
    except:
        raise TypeError('Something wrong with the input!')
    
    rs, thes = calc_polarcoords (center, dims )
    
    # anti distort
    for i in range(len(ns)):
        rs *= ncempy.algo.distortion.rad_dis(thes, dists[i*2+1], dists[i*2+2], ns[i]) 
    
    dis_xx = rs*np.cos(thes)+center[0]
    dis_yy = rs*np.sin(thes)+center[1]
    
    f = scipy.interpolate.RectBivariateSpline(dims[0][0], dims[1][0], img)

    return f.ev(dis_xx, dis_yy)
        
    
def calc_radialprofile( img, rs, rMax, dr, rsigma, mask=None):
    '''
    Calculate the radial profile using Gaussian density estimator.
    
    It is suggested to use an rMax such that all directions are still in the image, otherwise outer areas will contribute differently and lead to different signal-to-noise. A value of dr corresponding to 1/10 px and rsigma corresponding to 1 px are good parameters.
    
    input:
    - img       image to take radial profile of intensity from
    - rs        array containing the radial distance, same shape as img
    - rMax      maximum radial distance used for the radial profile
    - dr        stepsize for r-axis of radial distance
    - rsigma    sigma for Gaussian used as kernel density estimator
    - mask      binary image, 0 for pixels to exclude
    
    return:
    - r, intens
    '''
    
    # check input
    try:
        assert(isinstance(img, np.ndarray))
        assert(isinstance(rs, np.ndarray))
        assert(np.array_equal(img.shape, rs.shape))
        
        rMax = float(rMax)
        dr = float(dr)
        rsigma = float(rsigma)
        
        if not mask is None:
            assert(isinstance(mask, np.ndarray))
            assert(np.array_equal(img.shape, mask.shape))
        
    except:
        raise TypeError('Something wrong with input.')
    
    # process mask
    if mask is None:
        mask = np.ones(img.shape)
    mask = mask.astype('float64')
    mask[np.where( mask > 0 )] = 1
    mask[np.where( mask == 0 )] = np.nan
    
    
    # prepare radial axis for hist  
    rBins = np.arange(0, rMax, dr)
    
    radialIndex = np.round(rs/dr)+1
    rBinMax = len(rBins)
    sel = (radialIndex<=rBinMax)
    
    sel = np.logical_and(sel, np.logical_not(np.isnan(mask)) )
    
    signal = np.histogram(rs[sel], rBins, weights=img[sel])
    count = np.histogram(rs[sel], rBins, weights=np.ones(img[sel].shape))

    signal_sm = scipy.ndimage.filters.gaussian_filter1d(signal[0], rsigma/dr)
    count_sm = scipy.ndimage.filters.gaussian_filter1d(count[0], rsigma/dr)
    
    # masked regions lead to 0 in count_sm, divide produces nans, just ignore the warning    
    old_err_state = np.seterr(divide='ignore',invalid='ignore')    
    signal_sm = np.divide(signal_sm,count_sm)
    np.seterr(**old_err_state)    
        
    return rBins[:-1], signal_sm
    
    
def plot_radialprofile( r, intens, dims, show=False ):
    '''
    Plot radial profile.
    
    input:
    - r             r-axis of radial profile
    - intens        intensity-axis of radial profile
    - dims          dimensions of original image to read out units
    
    return:
    - plot          plot rendered to np.ndarray
    '''
    
    try:
        # check data
        assert(isinstance(r, np.ndarray))
        assert(isinstance(intens, np.ndarray))
        assert(np.array_equal(r.shape, intens.shape))
        
        # check if dims availabel
        assert(len(dims)>=1)
        assert(len(dims[0])==3)

    except:
        raise TypeError('Something wrong with the input!')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(r, intens, 'r-')
    
    # labels
    ax.set_xlabel('r /{}'.format(dims[0][2]))
    ax.set_ylabel('I /[a.u.]')
        
    if show:
        plt.show(block=False)
    
    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return plot


def residuals_fit( param, r, intens, funcs ):
    '''
    Residual function to fit radial profile with flexibility.
    
    The arguments for the single functions are automatically cut from the fit parameters.
    
    input:
    - param         fit parameters
    - r             r-axis
    - intens        intensity-axis
    - funcs         list of functions to include
        
    return:
    - residuals
    '''

    return intens - ncempy.algo.math.sum_functions(r, funcs, param)
    
    
def fit_radialprofile( r, intens, funcs, init_guess, maxfev=None ):
    '''
    Fit the radial profile.
    
    Convenience wrapper for fitting.
    
    input:
    - r             r-axis of radial profile
    - intens        intensity-axis of radial profile
    - funcs         list of functions
    - init_guess    initial guess for parameters of functions in funcs
    '''    
    
    try:
        # check data
        assert(isinstance(r, np.ndarray))
        assert(isinstance(intens, np.ndarray))
        assert(np.array_equal(r.shape, intens.shape))
        
        # funcs and params
        assert(len(funcs)>=1)
        for i in range(len(funcs)):
            assert(funcs[i] in ncempy.algo.math.lkp_funcs)
        
        init_guess = np.array(init_guess)
        init_guess = np.reshape(init_guess, sum(map(lambda x: ncempy.algo.math.lkp_funcs[x][1], funcs)))

    except:
        raise TypeError('Something wrong with the input!')
    
    if maxfev is None:
        maxfev = 1000
 
    popt, flag = scipy.optimize.leastsq( residuals_fit, init_guess, args=(r, intens, funcs), maxfev=maxfev)

    if flag not in [1,2,3,4]:
        print('WARNING: fitting of radial profile failed.')

    return popt


def plot_fit( r, intens, dims, funcs, param, show=False ):
    '''
    Plot the fit results to the radial profile.
    
    input:
    - r             r-axis of radial profile
    - intens        intensity-axis of radial profile
    - dims          dimensions of original image to read out units
    - funcs         list of functions
    - param         parameters for functions in funcs
    
    return:
    - plot          plot rendered to np.ndarray
    '''
    
    try:
        # check data
        assert(isinstance(r, np.ndarray))
        assert(isinstance(intens, np.ndarray))
        assert(np.array_equal(r.shape, intens.shape))
        
        # check if dims availabel
        assert(len(dims)>=1)
        assert(len(dims[0])==3)
        
        # funcs and params
        assert(len(funcs)>=1)
        for i in range(len(funcs)):
            assert(funcs[i] in ncempy.algo.math.lkp_funcs)
            
        param = np.array(param)
        param = np.reshape(param, sum(map(lambda x: ncempy.algo.math.lkp_funcs[x][1], funcs)))

    except:
        raise TypeError('Something wrong with the input!')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # plot radial profile
    ax.plot(r, intens, 'r-')
    
    # plot single
    n = 0
    for i in range(len(funcs)):
        ax.plot(r, ncempy.algo.math.lkp_funcs[funcs[i]][0]( r, param[n:n+ncempy.algo.math.lkp_funcs[funcs[i]][1]]), 'g-' )
        n += ncempy.algo.math.lkp_funcs[funcs[i]][1]
    # sum of functions
    sum_funcs = ncempy.algo.math.sum_functions( r, funcs, param )
    ax.plot(r, sum_funcs, 'b-')
    
    # labels
    ax.set_xlabel('r /{}'.format(dims[0][2]))
    ax.set_ylabel('I /[a.u.]')
        
    if show:
        plt.show(block=False)
    
    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return plot
    
    
def run_singleImage( img, dims, settings, show=False ):
    '''
    Evaluate a single ring diffraction pattern with given settings.
    
    input:
    - img           image as np.ndarray
    - dims          corresponding dim vectors
    - settings      dict of settings necessary for the evaluation
        - lmax_r        (int)               radius used for local maximum detection                 [px]
        - lmax_thresh   (int)               threshold for local maximum detection                   [intensity]
        - lmax_cinit    (int, int)          initial guess for center                                [px]
        - lmax_range    (float, float)      r range to cinit used to filter local maxima            [dims]
        - ns            (int, ...)          distortion orders to correct                            []
        - fit_rrange    (float, float)      r range used to fit radial profile                      [dims]
        - back_xs       (float, ...)        fixpoints for subtracting background                    [dims]
        - back_xswidth  (float)             range around fixpoints to take into account             [dims]
        - back_init     (float)             initial guess for subtracting background
        - fit_funcs     (str, ...)          list of functions to model radial profile               []
        - fit_init      (float, ...)        initial guess for fitting                               []

        optional (set to None to use defaults)
        - plt_imgminmax (float, float)      relative range to plot the img                          []
        - rad_rmax      (float)             max r to be used in radial profile                      [dims]
        - rad_dr        (float)             stepsize for r-axis in radial profile                   [dims]
        - rad_sigma     (float)             sigma for Gaussian used as kernel density estimator     [dims]
        - mask          (np.ndarray)        binary image as img, 0 for pixels to exclude            []
        - fit_maxfev    (int)               maxfev forwarded to scipy optimize                      []
        
    return:
    - res
    '''
    
    try:
        assert(isinstance(img, np.ndarray))
        
        # check if dims availabel
        assert(len(dims)>=1)
        assert(len(dims[0])==3)
        
        assert(type(settings) is dict)
        
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
            mysettings['plt_imgminmax'] = (0.,1.)
            
        plot = ncempy.algo.local_max.plot_points(img, points, vminmax=mysettings['plt_imgminmax'], dims=dims, invert=True, show=show)
    
    # optimize center
    center = ncempy.algo.distortion.optimize_center(points, center_init, verbose=show)
        
    # fit distortions
    points_plr = ncempy.algo.distortion.points_topolar(points, center)
    dists = ncempy.algo.distortion.optimize_distortion(points_plr, mysettings['ns'])
    if show:
        plot = ncempy.algo.distortion.plot_distpolar(points_plr, dims, dists, mysettings['ns'], show=show)
    
    # calc coordinates in optimized system
    rs, thes = ncempy.algo.radial_profile.calc_polarcoords( center, dims, mysettings['ns'], dists )
    
    if settings['rad_rmax'] is None:
        mysettings['rad_rmax'] = np.abs(dims[0][0][0]-dims[0][0][1])*np.min(img.shape)/2.0
    if settings['rad_dr'] is None:
        mysettings['rad_dr'] = np.abs(dims[0][0][0]-dims[0][0][1])/10.
    if settings['rad_sigma'] is None:
        mysettings['rad_sigma'] = np.abs(dims[0][0][0]-dims[0][0][1])
        
    # extract radial profile
    R, I = ncempy.algo.radial_profile.calc_radialprofile( img, rs, mysettings['rad_rmax'], mysettings['rad_dr'], mysettings['rad_sigma'], mask=mysettings['mask'] )
    
    # cut radial profile
    sel = (R>=mysettings['fit_rrange'][0])*(R<=mysettings['fit_rrange'][1])
    I = I[sel]
    R = R[sel]
    
    rawRI = np.copy(np.array([R,I]).transpose())
    
    # subtract a power law background fitted to specific points
    fit_R = np.array([])
    fit_I = np.array([])
    for xpoint in mysettings['back_xs']:
        ix = np.where( np.abs(R-xpoint) < mysettings['back_xswidth'] )
        fit_R = np.append(fit_R, R[ix])
        fit_I = np.append(fit_I, I[ix])
            
    funcs_back = [ 'const', 'powlaw']
    res_back = ncempy.algo.radial_profile.fit_radialprofile( fit_R, fit_I, funcs_back, mysettings['back_init'], maxfev=1000 )
    if show:
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, funcs_back, res_back, show=show )
        
    I = I - ncempy.algo.math.sum_functions(R, funcs_back, res_back)    
    
    
    # fit
    res = ncempy.algo.radial_profile.fit_radialprofile( R, I, mysettings['fit_funcs'], mysettings['fit_init'], maxfev=mysettings['fit_maxfev'])
    
    if show:
        plot = ncempy.algo.radial_profile.plot_fit( R, I, dims, mysettings['fit_funcs'], res, show=show )
        
    return np.array([R,I]).transpose(), res, center, dists, rawRI, res_back, mysettings
    
