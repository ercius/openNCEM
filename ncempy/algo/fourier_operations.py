import numpy as np


def shearImage(im, dim, shfactor):
    """ Exact shear of an image using the frequency space of the FFT
    of the image

    Parameters
    ----------
        im: ndarray, 2D
            The image to shear.
        dim: int
            The axis to shear along
        shfactor: float
            The amount of shear.
    Returns
    -------
        : ndarray
            The sheared image.
    """
    fftF = np.fft.fftshift(np.fft.fftfreq(im.shape[dim], d=1 / im.shape[dim]))
    # xx, yy = np.meshgrid(np.linspace(0,im.shape[0]-1,im.shape[0]),
    #                     np.linspace(0,im.shape[1]-1,im.shape[1]))
    xx, yy = np.meshgrid(fftF, fftF)
    f = np.fft.fftshift(np.fft.fftn(im, axes=(dim,)), axes=(dim,))
    ff = f * np.exp(-2.0j * np.pi / im.shape[dim] * xx * yy * shfactor)

    sh = np.fft.ifftn(np.fft.ifftshift(ff, axes=(dim,)), axes=(dim,))
    return sh


def shiftImage(im, shift):
    """ Exact shear of an image using the frequency space of the FFT
    of the image

    Parameters
    ----------
        im: ndarray, 2D
            The image to shear.
        shift: tuple of floats
            The number of pixels to shift in each direction
    Returns
    -------
        : ndarray
            The shifted image.
    """
    fftF0 = np.fft.fftshift(np.fft.fftfreq(im.shape[0], d=1 / im.shape[0]))
    fftF1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d=1 / im.shape[1]))
    xx, yy = np.meshgrid(fftF0, fftF1)
    f = np.fft.fftshift(np.fft.fftn(im))
    ff = f * np.exp(-2.0j * np.pi / im.shape[0] * xx * shift[0] +
                    -2.0j * np.pi / im.shape[1] * yy * shift[1])

    sh = np.fft.ifftn(np.fft.ifftshift(ff))
    return sh


def rotateImage(im, theta):
    """ Use three shears in Fourier space to exactly (and reversibly) rotate
    an image

    alpha = gamma = -tan(theta/2)
    beta = sin(theta)

    Need three shears: shearX, shearY, shearX
    three shears = [1+alpha*beta, alpha+gamma+alpha*beta*gamma; beta, 1 + beta*gamma]

    im1 = shearImage(im, 0, alpha)
    im2 = shearImage(im1, 1, beta)
    im3 = shearImage(im2,0, gamma)

    """
    print('Not implemented yet')
