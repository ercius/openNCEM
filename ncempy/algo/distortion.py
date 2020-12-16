"""
Module to handle distortions in diffraction patterns.
"""

import numpy as np
import scipy.optimize


def filter_ring(points, center, rminmax):
    """Filter points to be in a certain radial distance range from center.

    Parameters
    ----------
        points : np.ndarray
            Candidate points.
        center : np.ndarray or tuple
            Center position.
        rminmax : tuple
            Tuple of min and max radial distance.

    Returns
    -------
        : np.ndarray
            List of filtered points, two column array.

    """
    
    try:
        # points have to be 2D array with 2 columns
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        assert(len(points.shape) == 2)    
        
        # center can either be tuple or np.array    
        center = np.array(center)
        center = np.reshape(center, 2)
        
        rminmax = np.array(rminmax)
        rminmax = np.reshape(rminmax, 2)
   
    except:
        raise TypeError('Something wrong with the input!')
    
    # calculate radii
    rs = np.sqrt( np.square(points[:,0]-center[0]) + np.square(points[:,1]-center[1]) )
    
    # filter by given limits
    sel = (rs>=rminmax[0])*(rs<=rminmax[1])
    
    if sel.any():
        return points[sel]
    else:
        return None
    
    
def points_topolar(points, center):
    """Convert points to polar coordinate system.

    Can be either in pixel or real dim, but should be the same for points and center.

    Parameters
    ----------
        points : np.ndarray
            Positions as two column array.
        center : np.ndarray or tuple
            Origin of the polar coordinate system.

    Returns
    -------
        : np.ndarray
            Positions in polar coordinate system as two column array (r, theta).

    """
    
    try:
        # points have to be 2D array with 2 columns
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        assert(len(points.shape) == 2)    
        
        # center can either be tuple or np.array    
        center = np.array(center)
        center = np.reshape(center, 2)
   
    except:
        raise TypeError('Something wrong with the input!')
    
    # calculate radii
    rs = np.sqrt( np.square(points[:,0]-center[0]) + np.square(points[:,1]-center[1]) )
    # calculate angle
    thes = np.arctan2(points[:,1]-center[1], points[:,0]-center[0])
    
    return np.array( [rs, thes] ).transpose()


def residuals_center( param, data):
    """Residual function for minimizing the deviations from the mean radial distance.

    Parameters
    ----------
        param : np.ndarray
            The center to optimize.
        data : np.ndarray
            The points in x,y coordinates of the original image.

    Returns
    -------
        : np.ndarray
           Residuals.

    """
    
    # manually calculating the radii, as we do not need the thetas
    rs = np.sqrt( np.square(data[:,0]-param[0]) + np.square(data[:,1]-param[1]) )
    
    return rs-np.mean(rs)
    
    
def optimize_center(points, center, maxfev=1000, verbose=None):
    """Optimize the center by minimizing the sum of square deviations from the mean radial distance.

    Parameters
    ----------
        points : np.ndarray
            The points to which the optimization is done (x,y coords in org image).
        center : np.ndarray or tuple
            Initial center guess.
        maxfev : int
            Max number of iterations forwarded to scipy.optimize.leastsq().
        verbose : bool
            Set to get verbose output.

    Returns
    -------
        : np.ndarray
            The optimized center.

    """

    try:
        # points have to be 2D array with 2 columns
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        assert(len(points.shape) == 2)    
        
        # center can either be tuple or np.array    
        center = np.array(center)
        center = np.reshape(center, 2)
   
    except:
        raise TypeError('Something wrong with the input!')

    # run the optimization
    popt, flag = scipy.optimize.leastsq(residuals_center, center, args=points, maxfev=maxfev)
    
    if flag not in [1, 2, 3, 4]:
        print('WARNING: center optimization failed.')

    if verbose:
        print('optimized center: ({}, {})'.format(center[0], center[1]))

    return popt
    

def rad_dis(theta, alpha, beta, order=2):
    """Radial distortion due to ellipticity or higher order distortion.

    Relative distortion, to be multiplied with radial distance.

    Parameters
    ----------
        theta : np.ndarray
            Angles at which to evaluate. Must be float.
        alpha : float
            Orientation of major axis.
        beta : float
            Strength of distortion (beta = (1-r_min/r_max)/(1+r_min/r_max).
        order : int
            Order of distortion.

    Returns
    -------
        : np.ndarray
            Distortion factor.
    """
    
    return (1.-np.square(beta))/np.sqrt(1.+np.square(beta)-2.*beta*np.cos(order*(theta+alpha)))

    
def residuals_dis(param, points, ns):
    """Residual function for distortions.

    Parameters
    ----------
        param : np.ndarray
            Parameters for distortion.
        points : np.ndarray
            Points to fit to.
        ns : tuple
            List of orders to account for.

    Returns
    -------
        : np.ndarray
           Residuals.
    """

    est = param[0]*np.ones(points[:, 1].shape)
    for i in range(len(ns)):
        est *=rad_dis(points[:, 1], param[i*2+1], param[i*2+2], ns[i])
            
    return points[:, 0] - est
    
    
def optimize_distortion(points, ns, maxfev=1000, verbose=False):
    """Optimize distortions.

    The orders in the list ns are first fitted subsequently and the result is refined in a final fit simultaneously fitting all orders.

    Parameters
    ----------
        points : np.ndarray
            Points to optimize to (in polar coords).
        ns : tuple
            List of orders to correct for.
        maxfev : int
            Max number of iterations forwarded to scipy.optimize.leastsq().
        verbose : bool
            Set for verbose output.

    Returns
    -------
        : np.ndarray
            Optimized parameters according to ns.
    """
    
    try:
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        # check points to be sufficient for fitting
        assert(points.shape[0] >= 3)
        
        # check orders
        assert(len(ns)>=1)
    except:
        raise TypeError('Something wrong with the input!')

    # init guess for full fit
    init_guess = np.ones(len(ns)*2+1)
    init_guess[0] = np.mean(points[:,0])
    
    # make a temporary copy
    points_tmp = np.copy(points)
    
    if verbose:
        print('correction for {} order distortions.'.format(ns))
        print('starting with subsequent fitting:')
    
    # subsequently fit the orders
    for i in range(len(ns)):
        # optimize order to points_tmp
        popt, flag = scipy.optimize.leastsq(residuals_dis, np.array((init_guess[0], 0.1, 0.1)),
                                                     args=(points_tmp, (ns[i],)), maxfev=maxfev)
        
        if flag not in [1, 2, 3, 4]:
            print('WARNING: optimization of distortions failed.')
        
        # information
        if verbose:
            print('fitted order {}: R={} alpha={} beta={}'.format(ns[i], popt[0], popt[1], popt[2]))
        
        # save for full fit
        init_guess[i*2+1] = popt[1]
        init_guess[i*2+2] = popt[2]
        
        # do correction
        points_tmp[:, 0] /= rad_dis(points_tmp[:, 1], popt[1], popt[2], ns[i])
    
    # full fit    
    if verbose:
        print('starting the full fit:')    
    
    popt, flag = scipy.optimize.leastsq(residuals_dis, init_guess, args=(points, ns), maxfev=maxfev)
    
    if flag not in [1, 2, 3, 4]:
        print('WARNING: optimization of distortions failed.')
    
    if verbose:
        print('fitted to: R={}'.format(popt[0]))
        for i in range(len(ns)):
            print('.. order={}, alpha={}, beta={}'.format(ns[i], popt[i*2+1], popt[i*2+2]))

    return popt
