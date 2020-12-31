import numpy as np


def moments(im, order=3):
    """ Calculate the image raw moments.

    Parameters
    ----------
        im : ndarray
            A 2D array of the image to calculate the moments for

        order : int, default=3
            The order to calculate the moments to.

    Returns
    -------
        : ndarray
            An ndarray of (order+1, order+1)
    """
    order = int(order) + 1
    XX, YY = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    mm = np.zeros((order, order))
    for ii in range(order):
        XXii = XX ** ii
        for jj in range(order):
            mm[ii, jj] = np.sum(XXii * YY ** jj * im.astype(np.float64))
    return mm


def moments_central(im, cent=None, order=3):
    """ Calculate the image central moments.

    Parameters
    ----------
        im : ndarray
            A 2D array of the image to calculate the central moments for

        cent : tuple, 2 values, optional
            The centroid calculated from the image raw moments. If it
            is not supplied then it is calculated by moments()

        order : int, default = 3
            The order to calculate the orders to.

    Returns
    -------
        : ndarray
            An ndarray of size (order + 1, order +1)

    """
    if not cent:
        m = moments(im, order)
        cent = centroid(m)

    order = int(order) + 1

    XX, YY = np.mgrid[-cent[0]:im.shape[0] - cent[0],
                      - cent[1]:im.shape[1] - cent[1]]
    XX = XX.astype('float64')
    YY = YY.astype('float64')

    mc = np.zeros((order, order))
    for ii in range(order):
        XXii = XX ** ii
        for jj in range(order):
            mc[ii, jj] = np.sum(XXii * YY ** jj * im.astype(np.float64))
    return mc


def centroid(m):
    """ Find the centroid based on the raw moments.

    Parameters
    ----------
        m : ndarray
            The output from moments()

    Returns
    -------
        : tuple, 2 elements (X,Y)
            The position of the centroid as a 2-tuple
    """
    return m[1, 0] / m[0, 0], m[0, 1] / m[0, 0]


def moment_angle(mc):
    """ Calculate the orientation of the object from its central moments

    Parameters
    ----------
        mc : ndarray
            The array of the central moments. Output from moments_central

    Returns
    -------
        : float
            The angle in radians.

    """
    mu_20 = mc[2, 0] / mc[0, 0]  # - centroid[0]**2
    mu_02 = mc[0, 2] / mc[0, 0]  # - centroid[1]**2
    mu_11 = mc[1, 1] / mc[0, 0]  # - centroid[0]*centroid[1]

    th = 0.5 * np.arctan2(2 * mu_11, mu_20 - mu_02)
    return th
