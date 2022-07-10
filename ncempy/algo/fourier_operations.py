import numpy as np


def shearImage(im, dim, shear_factor):
    """ Exact shear of an image using the frequency space of the FFT
    of the image.

    Currently only works for square images.

    Parameters
    ----------
        im: ndarray, 2D
            The image to shear.
        dim: int
            The axis to shear along
        shear_factor: float
            The amount of shear.

    Returns
    -------
        : ndarray
            The sheared image.
    """

    if im.shape[0] != im.shape[1]:
        print('Input image must be square.')
        return

    fftF = np.fft.fftshift(np.fft.fftfreq(im.shape[dim], d=1 / im.shape[dim]))
    yy, xx = np.meshgrid(fftF, fftF, indexing='ij')
    f = np.fft.fftshift(np.fft.fftn(im, axes=(dim,)), axes=(dim,))
    ff = f * np.exp(-2.0j * np.pi / im.shape[dim] * xx * yy * shear_factor)

    sh = np.fft.ifftn(np.fft.ifftshift(ff, axes=(dim,)), axes=(dim,))
    sh = np.real(sh)
    return sh


def shiftImage(im, shift):
    """ Exact shear of an image using the frequency space of the FFT
    of the image.

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
    xx, yy = np.meshgrid(fftF0, fftF1, indexing='ij')
    f = np.fft.fftshift(np.fft.fftn(im))
    ff = f * np.exp(-2.0j * np.pi / im.shape[0] * xx * shift[0] +
                    -2.0j * np.pi / im.shape[1] * yy * shift[1])

    sh_im = np.fft.ifftn(np.fft.ifftshift(ff))
    sh_im = np.real(sh_im)
    return sh_im


def rotateImage(im, theta, pad=False):
    """ Use three shears in Fourier space to exactly (and reversibly) rotate
    an image.

    Currently only works for square images.

    Parameters
    ----------
        im: ndarray, 2D
            The image to rotate

        theta: float
            The angle to rotate by in radians

        pad: bool
            Add padding to the image before rotating

    """

    if im.shape[0] != im.shape[1]:
        print('Input image must be square.')
        return

    alpha = gamma = - np.tan(theta/2.)
    beta = np.sin(theta)

    if pad:
        p = im.shape[0] * np.sin(theta)
        im = np.pad(im, int(p//2+1))

    im1 = shearImage(im,  0, alpha)
    im2 = shearImage(im1, 1, beta)
    im3 = shearImage(im2, 0, gamma)
    im3 = np.real(im3)
    return im3


def bandpass_filter(im, inner, outer, sigma):
    """ Apply a bandpass filter with smoothed edges in Fourier space.

    Parameters
    ----------
    im: ndarray, 2D
        The image to apply the bandpass filter to.
    inner: float
        The lowest spatial frequency allowed in terms of the Nyquist frequency. Should be between 0 and 0.5.
    outer: float
        The highest spatial frequency allowed in terms of the Nyquist frequency. Should be between 0 and 0.5.
    sigma: float
        The standard deviation of the Gaussian smoothing to remove hard edges.

    Returns
    -------
    : ndarray, 2D
        The image with the bandpass filter applied.

    """
    # Create 2D fft frequency array which goes from 0 to 0.5
    num0, num1 = im.shape
    nf0 = np.fft.fftshift(np.fft.fftfreq(num0))
    nf1 = np.fft.fftshift(np.fft.fftfreq(num1))
    XX, YY = np.meshgrid(nf0, nf1, indexing='ij')

    # Create a mask that is 1 between the two desired frequencies.
    mask = np.sqrt(XX ** 2 + YY ** 2)
    mask = (mask < outer) * (mask > inner)  # bool data type

    # Gaussian smooth the edges to remove edge artifacts
    # Gaussian peak is at the [num/2,num/2] pixel; important for phase
    linG0 = np.linspace(-(num0 / 2), num0 / 2 - 1, num0)
    linG1 = np.linspace(-(num1 / 2), num1 / 2 - 1, num1)
    gXX, gYY = np.meshgrid(linG0, linG1, indexing='ij')
    gg = np.exp(-(gXX ** 2 / (2 * sigma ** 2) + gYY ** 2 / (2 * sigma ** 2))) / (
                2 * np.pi * sigma * sigma)  # 2D gaussian with unit area

    # Apply the smoothing to the mask using convolution
    gg_fft = np.fft.rfft2(gg)
    mask_fft = np.fft.rfft2(mask)
    cc = gg_fft * mask_fft
    mask_smooth = np.fft.irfft2(cc)

    # Apply the mask to the image in Fourier space
    gg_filt_fft = np.fft.fft2(im) * mask_smooth
    gg_filt = np.real(np.fft.ifft2(gg_filt_fft))

    return gg_filt
