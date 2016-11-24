'''
Module to find local maxima in an image.
'''

import numpy as np
import scipy.ndimage.filters

import matplotlib.pyplot as plt


def local_max(img, r, thresh):
    '''Find local maxima from comparing dilated and eroded images.
    
    Calculates images with maximum and minimum within given radius. If the difference is larger than the threshold, the original pixel position with max value is detected as local maximum.
    
    Parameters:
        img (np.ndarray):    Input image.
        r (int):    Radius for locality.
        thresh (int/float):    Intensity difference threshold.
    
    Returns:
        (np.ndarray):    Array of points.
        
    '''
    
    try:
        r = int(r)
        thresh = int(thresh)
        assert(isinstance(img, np.ndarray))
    except:
        raise TypeError('Bad input!')
    
    
    # prepare circular kernel
    y,x = np.ogrid[-r:r+1, -r:r+1]
    kernel = x**2 + y**2 <= r**2
    
    # calculate max and min images
    img_dil = scipy.ndimage.filters.maximum_filter(img, footprint=kernel)
    img_ero = scipy.ndimage.filters.minimum_filter(img, footprint=kernel)
    
    # get selection of local maxima
    sel = (img==img_dil)*(img-img_ero > thresh)
    
    if sel.any():  
        # retrieve and return points
        points = np.argwhere(sel)
        #points = np.roll(points, 1, axis=1)
        
        return points
    else:
        # otherwise return None to avoid having an empty list
        return None
        
    
def plot_points(img, points, vminmax=(0,1), dims=None, invert=False, show=False):
    '''Plot the detected points on the input image for checking.
    
    Parameters:
        img (np.ndarray):    Image.
        points (np.ndarray):    Array containing the points.
        vminmax (tuple):    Tuple of two values for relative lower and upper cut off to display image.
        dims (tuple):    Tuple of dims to plot in dimensions.
        invert (bool):    Set to invert the image.
        show (bool):    Set to directly show the plot interactively.
    
    Returns:
        (np.ndarray):    Image of the plot.
        
    '''
    
    try:
        assert(isinstance(img, np.ndarray))
        
        assert(isinstance(points, np.ndarray))
        assert(points.shape[1] == 2)
        assert(len(points.shape) == 2)
    except:
        raise TypeError('Something wrong with the input!')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if invert:
        cmap="Greys"
    else:
        cmap="gray"
    
    if dims:
        ax.imshow(img, cmap=cmap, vmin=np.min(img)+vminmax[0]*(np.max(img)-np.min(img)), vmax=np.min(img)+vminmax[1]*(np.max(img)-np.min(img)), extent=(np.min(dims[0][0]), np.max(dims[0][0]), np.max(dims[1][0]), np.min(dims[1][0])) )
        ax.set_xlabel('{} {}'.format(dims[0][1], dims[0][2]))
        ax.set_ylabel('{} {}'.format(dims[1][1], dims[1][2]))
        ax.set_xlim((np.min(dims[0][0]), np.max(dims[0][0])))
        ax.set_ylim((np.max(dims[1][0]), np.min(dims[1][0])))
    else:
        ax.imshow(img, cmap=cmap, vmin=np.min(img)+vminmax[0]*(np.max(img)-np.min(img)), vmax=np.min(img)+vminmax[1]*(np.max(img)-np.min(img)) )
        ax.set_xlim((0,img.shape[1]-1))
        ax.set_ylim((img.shape[0]-1,0))
    
    ax.scatter(points[:,1], points[:,0], color='r', marker='o', facecolors='none')
    
    if show:
        try:
            plt.show(block=False)
        except:
            pass
    
    fig.canvas.draw()
    
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return plot
    
    
def points_todim(points, dims):
    '''Convert points from px coordinates to real dim.
    
    Points are expected to be array indices for the first two dimensions in dims.
    
    Parameters:
        points (np.ndarray):    Points to convert.
        dims (tuple):    Tuple of dimensions.
        
    Returns:
        (np.ndarray):   Converted points.
        
    '''
    
    try:
        # try to convert input to np.ndarray with 2 columns (necessary if only one entry provided)
        points = np.reshape(np.array(points), (-1,2))
        # check if enough dims availabel
        assert(len(dims)>=2)
        assert(len(dims[0])==3)
    except:
        raise TypeError('Something wrong with the input!')
    
    # do the conversion by looking up thing in dimension vectors
    points_d = np.array( [ dims[0][0][points[:,0]], dims[1][0][points[:,1]] ] ).transpose()
    
    return points_d
