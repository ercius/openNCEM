""" A module with Gaussian functions for creating Gaussian distributions
in 1D, 2D and 3D. There are also functions with the same name and _FIT
which can be used to fit an intensity distribution in a ndarray with the
corresponding Gaussian distribution using scipy's curve_fit.

A 1D Lorentz and a 1D gaussLorentz function is also included for fitting
zero loss peaks in EELS spectra.

These functions were written assuming you use np.meshgrid() with indexing = 'xy'. If
you use 'ij' indexing then be careful about how you pass in the x and y coordinates.

"""

import numpy as np


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
            A vector of size (N,) of the Gaussian distribution evaluated
            at the input x values.
            
    Note
    ----
        Calculate the half width at half maximum (HWHM) as
        >> HWHM = sqrt(2*log(2))*stDev ~ 1.18*stDev
        or if x goes from -1 to 1
        >> 0.5*size(x)*stDev
        
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
    return w ** 2 / ((x - x0) ** 2 + w ** 2)


def gaussLorentz1D(x, x0, w):
    """ A Gaussian-Lorentzian function in one dimension.

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
        lg: ndarray
            A vector of size (N,) of the Lorentzian Gaussian distribution
            evaluated at the input x values.

    """
    gSigma = w / (2 * np.sqrt(2 * np.log(2)))  # change the width to sigma
    return gauss1D(x, x0, gSigma) + lorentz1D(x, x0, w)


def gauss2D(x, y, x0, y0, sigma_x, sigma_y):
    """ Calculates the value of a Gaussian at a 2D set of points for the given
    standard deviations with maximum normalized to 1.
    The Gaussian axes are assumed to be 90 degrees from each other.
    
    Parameters
    ----------
        x, y : ndarray
            A 2D array of size (M,N) of points for the x values. 
            Use x, y = np.meshgrid(range(M),range(N),indexing='xy').
        x0 : float
            The center of the Gaussian along the x direction.
        y0 : float
            The center of the Gaussian along the y direction.
        sigma_x : float
            Standard deviation of the Gaussian along x.
        sigma_y : float
            Standard deviation of the Gaussian along y.
            
    Returns
    -------
        g : ndarray
            A ndarray of size (N, M) of the Gaussian distribution evaluated
            at the input x values.
            
    Note
    ----
        The Gaussian is normalized such that the peak == 1. To normalize the
        integral divide by 2 * np.pi * sigma_x * sigma_y
        
    """
    x0 = float(x0)
    y0 = float(y0)
    g2 = np.exp(-((x - x0) ** 2 / (2 * sigma_x ** 2) + (y - y0) ** 2 / (2 * sigma_y ** 2)))
    g2_norm = g2 / np.max(g2.flatten())
    return g2_norm


def gauss2D_FIT(xy, x0, y0, sigma_x, sigma_y):
    """ Version of gauss2D used for fitting (1 x_data input (xy)
    and flattened output). Returns the value of a gaussian at a 2D set of points for the given
    standard deviation with maximum normalized to 1.
    The Gaussian axes are assumed to be 90 degrees from each other.
    
    Parameters
    ----------
        xy: tuple
            A (N,M) shaped tuple containing the vectors of points for the
            evaluation points.
        x0: float
            The x center of the Gaussian.
        y0: float
            The y center of the Gaussian.
        sigma_x: float
            The standard deviation of the Gaussian along x.
        sigma_y: float
            The standard deviation of the Gaussian along y.
    
    Returns
    -------
        g2_norm: ndarray
            The Gaussian distribution with maximum value normalized to 1. The
            2D ndarray is reshaped into a (N*M,) array for use in fitting
            functions in numpy and scipy.
        
    """
    x0 = float(x0)
    y0 = float(y0)
    x = xy[0]
    y = xy[1]
    g2 = np.exp(-((x - x0) ** 2 / (2 * sigma_x ** 2) + (y - y0) ** 2 / (2 * sigma_y ** 2)))
    g2_norm = g2 / np.max(g2.flatten())
    return g2_norm.reshape(-1)


def gauss2D_theta(x, y, x0, y0, sigma_x, sigma_y, theta):
    """ Returns the value of a generalized Gaussian at a 2D set of points for
    the given standard deviation with maximum normalized to 1.
    The Gaussian axes (assumed to be 90 degrees from each other) can be
    oriented at different angles to the output array axes using theta.
    
    Parameters
    ----------
        x: ndarray or list
            Evaluation points along x of size (N,)
        y: ndarray
            Evaluation points along y of size (M,)
        x0: float
            Center of the maximum along x.
        y0: float
            Center of the maximum along y.
        sigma_x: float
            Standard deviation along x.
        sigma_y: float
            Standard deviation along y.
        theta: float
            The rotation of the Gaussian principle axes from the array
            horizontal and vertical. In radians.
            
    Returns
    -------
        g2_norm: ndarray
            A (N, M) sized 2D ndarray with a Gaussian distribution rotated.
            
    """

    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g2 = np.exp(- (a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) ** 2)))
    g2_norm = g2 / np.max(g2.flatten())
    return g2_norm  # return a 2D array


def gauss2D_theta_FIT(xy, x0, y0, sigma_x, sigma_y, theta):
    """ Version of gauss2D_theta used for fitting (1 x_data input and
    flattened output). Returns the value of a gaussian at a 2D set of points for the given standard deviation with
    maximum normalized to 1. The Gaussian axes can be oriented at different angles (theta) in radians.
    
    Parameters
    ----------
        xy: tuple
            Evaluation points along x and y of shape (N,M)
        x0: float
            Center of the maximum along x.
        y0: float
            Center of the maximum along y.
        sigma_x: float
            Standard deviation along x.
        sigma_y: float
            Standard deviation along y.
        theta: float
            The rotation of the Gaussian principle axes from the array
            horizontal and vertical. In radians.
            
    Returns
    -------
        g2_norm: ndarray
            A (N, M) sized 2D ndarray with a Gaussian distribution rotated.

    """
    x0 = float(x0)
    y0 = float(y0)
    x = xy[0]
    y = xy[1]
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g2 = np.exp(- (a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) ** 2)))
    g2_norm = g2 / np.max(g2.flatten())
    return g2_norm.ravel()  # return a 1D vector


def gauss2D_poly_FIT(xy, x0, y0, A, B, C):
    """ Returns the flattened values of an elliptical gaussian at a 2D set of
    points for the given polynomial pre-factors with maximum normalized to 1.
    The Gaussian axes are assumed to be 90 degrees from each other.
    
    The matrix looks like [[A,B],[B,C]].
    See https://en.wikipedia.org/wiki/Gaussian_function
    
    Parameters
    ----------
        xy: tuple
            Evaluation points along x and y of shape (N,M)
        x0: float
            Center of the maximum along x.
        y0: float
            Center of the maximum along y.
        A: float
            The A pre-factor for the polynomial expansion in the exponent.
        B: float
            The B pre-factor for the polynomial expansion in the exponent.
        C: float
            The C pre-factor for the polynomial expansion in the exponent.
            
    Returns
    -------
        g2_norm: ndarray
            A (N, M) sized 2D ndarray with a Gaussaian distribution rotated.
    """

    x0 = float(x0)
    y0 = float(y0)
    x = xy[0]  # retrieve the array from the tuple
    y = xy[1]
    g2 = np.exp(-(A * (x - x0) ** 2 + 2 * B * (x - x0) * (y - y0) + C * (y - y0) ** 2))
    g2_norm = g2 / np.max(g2.flatten())
    return g2_norm.ravel()


def gauss3D(x, y, z, x0, y0, z0, sigma_x, sigma_y, sigma_z):
    """
    gauss3D(x,y,z,x0,y0,z0,sigma_x,sigma_y,sigma_z)
    Returns the value of a gaussian at a 2D set of points for the given
    standard deviations with maximum normalized to 1.
    The Gaussian axes are assumed to be 90 degrees from each other.

    Note
    -----
    Be careful about the indexing used in meshgrid and the order in which you pass the x, y, z variables in.

    Parameters
    ----------
        x, y, z: ndarray, from numpy.meshgrid
            2D arrays of points (from meshgrid)
     
        x0, y0, z0: float
            The x, y, z centers of the Gaussian
        
        sigma_x, sigma_y, sigma_z: float
            The standard deviations of the Gaussian.
    Returns
    -------
        g3_norm: ndarray
            A 3D ndarray
            
    """
    x0 = float(x0)
    y0 = float(y0)
    z0 = float(z0)
    g3 = np.exp(-((x - x0) ** 2 / (2 * sigma_x ** 2) + (y - y0) ** 2 / (2 * sigma_y ** 2) + (z - z0) ** 2 / (
                2 * sigma_z ** 2)))
    g3_norm = g3 / np.max(g3.flatten())
    return g3_norm


def gauss3D_FIT(xyz, x0, y0, z0, sigma_x, sigma_y, sigma_z):
    """
    gauss3D_FIT((x,y,z),x0,y0,z0,sigma_x,sigma_y,sigma_z)
    Returns the value of a gaussian at a 2D set of points for the given
    standard deviations with maximum normalized to 1.
    The Gaussian axes are assumed to be 90 degrees from each other.
     xyz - 
     x0, y0, z0 = the x, y, z centers of the Gaussian
     sigma_x, sigma_y, sigma_z = The std. deviations of the Gaussian.

    Note
    -----
    Be careful about the indexing used in meshgrid and the order in which you pass the x, y, z variables in.

    Parameters
    ----------
        xyz: tuple of ndarrays
            A tuple containing the 3D arrays of points (from meshgrid)
        x0, y0, z0: float
            The x, y, z centers of the Gaussian
        sigma_x, sigma_y, sigma_z: float
            The standard deviations of the Gaussian.
    Returns
    -------
        g3_norm: ndarray
            A flattened array for fitting.
    
    """

    x0 = float(x0)
    y0 = float(y0)
    z0 = float(z0)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    g3 = np.exp(-((x - x0) ** 2 / (2 * sigma_x ** 2) + (y - y0) ** 2 / (2 * sigma_y ** 2) + (z - z0) ** 2 / (
                2 * sigma_z ** 2)))
    g3_norm = g3 / np.max(g3.flatten())
    return g3_norm.ravel()


def gauss3D_poly(x, y, z, x0, y0, z0, A, B, C, D, E, F):
    """
    gauss3Dpoly_FIT((x,y,z),x0,y0,z0,A,B,C,D,E,F)
    Returns the value of a gaussian at a 2D set of points for the given
    standard deviations with maximum normalized to 1.
    The Gaussian axes are not locked to be 90 degrees.

    Parameters
    ----------
        x, y, z : ndarray, 3D
            3D arrays of points (from meshgrid)

        x0, y0, z0 : float
            The x, y, z centers of the Gaussian

        A, B, C, D, E, F : float
            The polynomial values for the fit

    Returns
    -------
        : ndarray, 3D
            The 3D Gaussian
    """
    x0 = float(x0)
    y0 = float(y0)
    z0 = float(z0)
    g3 = np.exp(-(A * (x - x0) ** 2 + B * (y - y0) ** 2 + C * (z - z0) ** 2 + 2 * D * (x - x0) * (y - y0) + 2 * E * (
                x - x0) * (z - z0) + 2 * F * (y - y0) * (z - z0)))
    g3_norm = g3 / np.max(g3.flatten())
    return g3_norm


def gauss3D_poly_FIT(xyz, x0, y0, z0, A, B, C, D, E, F):
    """
    gauss3Dpoly_FIT((x,y,z),x0,y0,z0,A,B,C,D,E,F)
    Returns the value of a gaussian at a 2D set of points for the given
    standard deviations with maximum normalized to 1.
    The Guassian axes are not locked to be 90 degrees.

    Parameters
    ----------
        xyz : tuple of ndarrays
            3D arrays of points (from meshgrid) combined in a tuple

        x0, y0, z0 : float
            The x, y, z centers of the Gaussian

        A, B, C, D, E, F : float
            The polynomial values for the fit

    Returns
    -------
        : ndarray, 3D
            The 3D Gaussian
    """
    x0 = float(x0)
    y0 = float(y0)
    z0 = float(z0)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    g3 = np.exp(-(A * (x - x0) ** 2 + B * (y - y0) ** 2 + C * (z - z0) ** 2 + 2 * D * (x - x0) * (y - y0) + 2 * E * (
                x - x0) * (z - z0) + 2 * F * (y - y0) * (z - z0)))
    g3_norm = g3 / np.max(g3.flatten())
    return g3_norm.ravel()


def gauss3DGEN_FIT(xyz, x0, y0, z0, sigma_x, sigma_y, sigma_z, Angle1, Angle2, Angle3, BG, Height):
    """
    gauss3DGEN_FIT((x,y,z),x0,y0,z0,sigma_x,sigma_y,sigma_z,Angle1,Angle2,Angle3,BG,Height)
    Returns the value of a gaussian at a 3D set of points for the given
    sub-pixel positions with standard deviations, 3D Eular rotation angles,
    constant Background value, and Gaussian peak height.

    adapted from code by Yongsoo Yang, yongsoo.ysyang@gmail.com

    Note
    ----
        This is a work in progress. Needs more testing.

    Parameters
    ----------
        xyz : tuple of 3 np.ndarray
            3D arrays of points (from meshgrid) combined in a tuple

        x0, y0, z0 : float
            The x, y, z centers of the Gaussian
        sigma_x,sigma_y,sigma_z: float
            standard deviations along x,y,z direction before 3D angular rotation
        Angle1, Angle2, Angle3 : float, degrees
            Tait-Bryan angles in ZYX convention for 3D rotation in degrees
        BG : float
            Background
        Height: float
            The peak height of the Gaussian function

    Returns
    -------
        : ndarray, 3D
            The 3D Gaussian
    """
    # 3D vectors for each sampled positions
    print('Warning: this function is using Fortran ordered array for use in tomviz. Need to test this')
    # Todo: Remove Fortran ordering
    v = np.array([xyz[0].reshape(-1, order='F') - x0,
                  xyz[1].reshape(-1, order='F') - y0,
                  xyz[2].reshape(-1, order='F') - z0])

    # rotation axes for Tait-Bryan angles
    vector1 = np.array([0, 0, 1])
    rotmat1 = MatrixQuaternionRot(vector1, Angle1)

    vector2 = np.array([0, 1, 0])
    rotmat2 = MatrixQuaternionRot(vector2, Angle2)

    vector3 = np.array([1, 0, 0])
    rotmat3 = MatrixQuaternionRot(vector3, Angle3)

    # full rotation matrix
    # Todo : Remove dependency on np.matrix()
    rotMAT = np.matrix(rotmat3) * np.matrix(rotmat2) * np.matrix(rotmat1)

    # 3x3 matrix for applying sigmas
    D = np.matrix(np.array([[1. / (2 * sigma_x ** 2), 0, 0, ],
                            [0, 1. / (2 * sigma_y ** 2), 0],
                            [0, 0, 1. / (2 * sigma_z ** 2)]]))

    # apply 3D rotation to the sigma matrix
    WidthMat = np.transpose(rotMAT) * D * rotMAT

    # calculate 3D Gaussian
    RHS_calc = WidthMat * np.matrix(v)
    Result = Height * np.exp(-1 * np.sum(v * RHS_calc.A, axis=0)) + BG

    return Result


def MatrixQuaternionRot(vector, theta):
    """
    Quaternion tp rotate a given theta angle around the given vector.

    adapted from code by Yongsoo Yang, yongsoo.ysyang@gmail.com

    Parameters
    ----------
        vector : ndarray, 3-element
            A non-zero 3-element numpy array representing rotation axis
        theta : float, degrees
            Rotation angle in degrees

    Returns
    -------

     Author: Yongsoo Yang, Dept. of Physics and Astronomy, UCLA
             yongsoo.ysyang@gmail.com
    """
    theta = theta * np.pi / 180
    vector = vector / np.float(np.sqrt(np.dot(vector, vector)))
    w = np.cos(theta / 2)
    x = -np.sin(theta / 2) * vector[0]
    y = -np.sin(theta / 2) * vector[1]
    z = -np.sin(theta / 2) * vector[2]
    RotM = np.array([[1. - 2 * y ** 2. - 2 * z ** 2, 2. * x * y + 2 * w * z, 2. * x * z - 2. * w * y],
                     [2. * x * y - 2. * w * z, 1. - 2. * x ** 2 - 2. * z ** 2, 2. * y * z + 2. * w * x],
                     [2 * x * z + 2 * w * y, 2 * y * z - 2. * w * x, 1 - 2. * x ** 2 - 2. * y ** 2]])

    return RotM
