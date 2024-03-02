"""
Module to find local maxima in an image.
"""

import numpy as np
import scipy.ndimage


def local_max(img, r, thresh):
    """Find local maxima from comparing dilated and eroded images.

    Calculates images with maximum and minimum within given radius. If the difference is larger than the threshold, the
    original pixel position with max value is detected as local maximum.

    Parameters
    ----------
        img : np.ndarray
            Input image.
        r : int
            Radius for locality.
        thresh : int or float
            Intensity difference threshold.

    Returns
    -------
        : np.ndarray
            Array of points.

    """
    
    try:
        r = int(r)
        thresh = int(thresh)
        assert(isinstance(img, np.ndarray))
    except:
        raise TypeError('Bad input!')

    # prepare circular kernel
    y, x = np.ogrid[-r:r+1, -r:r+1]
    kernel = x**2 + y**2 <= r**2
    
    # calculate max and min images
    img_dil = scipy.ndimage.maximum_filter(img, footprint=kernel)
    img_ero = scipy.ndimage.minimum_filter(img, footprint=kernel)
    
    # get selection of local maxima
    sel = (img == img_dil)*(img-img_ero > thresh)
    
    if sel.any():  
        # retrieve and return points
        points = np.argwhere(sel)

        return points
    else:
        # otherwise return None to avoid having an empty list
        return None


def points_todim(points, dims):
    """Convert points from px coordinates to real dim.

    Points are expected to be array indices for the first two dimensions in dims.

    Parameters
    ----------
        points : np.ndarray
            Points to convert.
        dims : tuple
            Tuple of dimensions.

    Returns
    -------
        : np.ndarray
           Converted points.

    """
    
    try:
        # try to convert input to np.ndarray with 2 columns (necessary if only one entry provided)
        points = np.reshape(np.array(points), (-1,2))
        # check if enough dims available
        assert(len(dims) >= 2)
        assert(len(dims[0]) == 3)
    except:
        raise TypeError('Something wrong with the input!')
    
    # do the conversion by looking up thing in dimension vectors
    points_d = np.array([dims[0][0][points[:, 0]], dims[1][0][points[:, 1]]]).transpose()
    
    return points_d
