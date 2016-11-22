'''
Module to handle distortions in diffraction patters.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize


def filter_ring(points, center, rminmax):
    '''
    Filter points to be in a certain radial distance range from center.
    
    input:
    - points        candidate points
    - center        center position
    - rminmax       tuple of min and max radial distance
    
    return:
    - points        list of filtered points
    '''
    
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
    '''
    Convert points to polar coordinate system.
    
    Can be either in pixel or real dim, but should be the same for points and center.
    
    input:
    - points        positions as two column array
    - center        origin of the polar coordinate system
    
    return:
    - points        positions in polar coordinate system as two column array (r, theta)
    '''
    
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
    
    
def plot_ringpolar(points, dims, show=False):
    '''
    Plot points in polar coordinate system.
    
    input:
    - points        positions in polar coords
    - dims          dimension information to plot labels
    '''

    try:
        # try to convert input to np.ndarray with 2 columns (necessary if only one entry provided)
        points = np.reshape(np.array(points), (-1,2))
        # check if enough dims availabel
        assert(len(dims)>=2)
        assert(len(dims[0])==3)
    except:
        raise TypeError('Something wrong with the input!')
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # mean value as line
    ax.axhline(np.mean(points[:,0]), ls='--', c='k')
    
    # points
    ax.plot(points[:,1], points[:,0], 'rx')
    
    # labels
    ax.set_xlabel('theta /[rad]')
    ax.set_xlim( (-np.pi, np.pi) )
    ax.set_ylabel('r /{}'.format(dims[0][2]))
    
    if show:
        plt.show(block=False)
    
    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return plot


def residuals_center( param, data):
    '''
    Residual function for minimizing the deviations from the mean radial distance.
    
    input:
    - param     the center to optimize
    - data      the points in x,y coordinates of the original image
    '''
    
    # manually calculating the radii, as we do not need the thetas
    rs = np.sqrt( np.square(data[:,0]-param[0]) + np.square(data[:,1]-param[1]) )
    
    return (rs-np.mean(rs))
    
    
def optimize_center(points, center, maxfev=1000, verbose=None):
    '''
    Optimize the center by minimizing the sum of square deviations from the mean radial distance.
    
    input:
    - points        the points to which the optimization is done (x,y coords in org image)
    - center        initial center guess
    - maxfev        max number of iterations forwarded to scipy.optimize.leastsq()
    
    return:
    - opt_center    the optimized center
    '''

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
    popt, flag = scipy.optimize.leastsq( residuals_center, center, args=(points), maxfev=maxfev)
    
    if flag not in [1,2,3,4]:
        print('WARNING: center optimization failed.')

    if verbose:
        print('optimized center: ({}, {})'.format(center[0], center[1]))

    return popt
    

def rad_dis( theta, alpha, beta, order=2 ):
    '''
    Radial distortion due to ellipticity or higher order distortion.
    
    Relative distortion, to be multiplied with radial distance.
    
    input:
    - theta     angles at which to evaluate
    - alpha     orientation of major axis
    - beta      strength of distortion (beta = (1-r_min/r_max)/(1+r_min/r_max)
    - order     order of distortion
    '''
    
    return (1.-np.square(beta))/np.sqrt(1.+np.square(beta)-2.*beta*np.cos(order*(theta+alpha)))

    
def residuals_dis(param, points, ns):
    '''
    Residual function for distortions.
    
    input:
    - param         parameters for distortion
    - points        points to compare to
    - ns            list of orders to account for
    
    return:
    - residuals
    '''

    est = param[0]*np.ones(points[:,1].shape)
    for i in range(len(ns)):
        est *=rad_dis( points[:,1], param[i*2+1], param[i*2+2], ns[i])
            
    return points[:,0] - est 
    
    
def optimize_distortion(points, ns, maxfev=1000, verbose=False):
    '''
    Optimize distortions.
    
    The orders in the list ns are first fitted subsequently and the result is refined in a final fit simultaneously fitting all orders.
    
    input:
    - points        points to optimize to (in polar coords)
    - ns            list of orders to correct for
    - maxfev        max number of iterations forwarded to scipy.optimize.leastsq()
    
    return:
    - popt          optimized parameters according to ns
    '''
    
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
        popt, flag = scipy.optimize.leastsq( residuals_dis, (init_guess[0], 0.1, 0.1), args=(points_tmp, (ns[i],)), maxfev=maxfev)
        
        if flag not in [1,2,3,4]:
            print('WARNING: optimization of distortions failed.')
        
        # information
        if verbose:
            print('fitted order {}: R={} alpha={} beta={}'.format(ns[i], popt[0], popt[1], popt[2]))
        
        # save for full fit
        init_guess[i*2+1] = popt[1]
        init_guess[i*2+2] = popt[2]
        
        # do correction
        points_tmp[:,0] /= rad_dis(points_tmp[:,1], popt[1], popt[2], ns[i])
    
    # full fit    
    if verbose:
        print('starting the full fit:')    
    
    popt, flag = scipy.optimize.leastsq( residuals_dis, init_guess, args=(points, ns), maxfev=maxfev)
    
    if flag not in [1,2,3,4]:
        print('WARNING: optimization of distortions failed.')
    
    if verbose:
        print('fitted to: R={}'.format(popt[0]))
        for i in range(len(ns)):
            print('.. order={}, alpha={}, beta={}'.format(ns[i], popt[i*2+1], popt[i*2+2]))

    return popt
    
    
def plot_distpolar(points, dims, dists, ns, show=False):
    '''
    Plot the results of distortion fitting in polar coordinates.
    
    input:
    - points        points used to fit to in polar coords
    - dims          dimensions, necessary to have unit information
    - dists         results of dist fitting, length according to ns
    - ns            list of used orders
    '''
    
    try:
        # check points
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        
        # check if enough dims availabel
        assert(len(dims)>=2)
        assert(len(dims[0])==3)
        
        # check orders
        assert(len(ns)>=1)
        
        # check dists
        assert(dists.shape[0] == len(ns)*2+1)
    except:
        raise TypeError('Something wrong with the input!')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # stuff from the single orders
    ax.axhline(dists[0], ls='--', c='k')
    xpl_ell = np.linspace(-np.pi, np.pi, 100)
    for i in range(len(ns)):
        plt.plot( xpl_ell, dists[0]*rad_dis(xpl_ell, dists[i*2+1], dists[i*2+2], ns[i]), 'm--')

    # points before
    ax.plot(points[:,1], points[:,0], 'rx')
    
    # sum of all distorts
    sum_dists = np.ones(xpl_ell.shape)*dists[0]
    for i in range(len(ns)):
        sum_dists *= rad_dis(xpl_ell, dists[i*2+1], dists[i*2+2], ns[i])
    plt.plot( xpl_ell, sum_dists, 'b-' )
    
    # points after
    points_corr = np.copy(points)
    for i in range(len(ns)):
        points_corr[:,0] /= rad_dis(points[:,1], dists[i*2+1], dists[i*2+2], ns[i])
    plt.plot( points_corr[:,1], points_corr[:,0], 'gx')
    
    # labels
    ax.set_xlabel('theta /[rad]')
    ax.set_xlim( (-np.pi, np.pi) )
    ax.set_ylabel('r /{}'.format(dims[0][2]))
    
    if show:
        plt.show(block=False)
    
    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return plot

