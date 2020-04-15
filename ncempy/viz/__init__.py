"""
A set of visualization functions based on matplotlib.pyplot which are generally
useful for S/TEM data visualization.

"""

import numpy as np
import matplotlib.pyplot as plt


def imsd(im, vmin=-2, vmax=2, **kwargs):
    """Show an array as an image with intensities compared to the standard deviation of the data. Other
    keyword args are passed to pyplot.imshow(). vmin and vmax are by default set to -2 and 2 respectively
    which are usually good values to set for S/TEM data.

        Parameters
        ----------
            im : ndarray
                The image to show.

        Keywords
        --------
            vmin, vmax : float, defulat = -2, 2
                The vmin and vmax values to pass to imshow.

        Returns
        -------
            : matplotlib.pyplt.Figure
                The handle to the created figure
    """
    fg, ax = plt.subplots(1, 1)
    im2 = im - im.mean()
    im3 = im2 / np.std(im2)
    ax.imshow(im3, vmin=vmin, vmax=vmax, **kwargs)
    return fg

def imfft(im, d=1.0, ax=None):
    """ Show a 2D FFT as a diffractogram with log scaling applied and zero frequency
    fftshifted tp the center. A new figure is created or an axis can be specified.

    The diffracotgram is calculated from the original intensities (I) as

    .. math::
        1 + 0.001 * I ^2

    Parameters
    ----------
        im: ndarray
            The 2D fft of the diffraction pattern to display as a diffractogram
        d: float, optional, default = 1.0
            The real space pixel size of the image used to get the FFT
        ax: pyplot axis, optional
            An axis to plot into.
    Returns
    -------
        : pyplot Figure
            The figure used to plot the diffractogram.

    Example
    -------
        This example shows how to display a 2D ndarray (image) as a
        diffractogram. The image has a real space pixel size of 0.1 nanometer.

        >>> imageFFT = np.fft.fft2(image)
        >>> impy.imfft(imageFFT, d = 0.1)

    """

    fftFreq0 = np.fft.fftshift(np.fft.fftfreq(im.shape[0], d))
    fftFreq1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d))
    if ax is None:
        fg, ax = plt.subplots(1, 1)
    ax.imshow(np.fft.fftshift(np.log(1 + .001 * np.abs(im) ** 2)),
              extent=(fftFreq0[0], fftFreq0[-1], fftFreq1[-1], fftFreq1[0]))
    return fg


def imrfft(im, d=1.0, ax=None):
    """Show a 2D rFFT (real FFT) as a diffractogram with log scaling applied
    and fftshifted along axis 0. See imfft for full details.

    Parameters
    ----------
        im: ndarray
            The 2D fft of the diffraction pattern to display as a diffractogram
        d: float, optional, default = 1.0
            The real space pixel size of the image used to get the FFT
        ax: pyplot axis, optional
            An axis to plot into.
    Returns
    -------
        : pyplot Figure
            The figure used to plot the diffractogram.
    """
    fftFreq1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d))
    fftFreq0 = np.fft.rfftfreq(im.shape[0], d)
    if ax is None:
        fg, ax = plt.subplots(1, 1)
    ax.imshow(np.fft.fftshift(np.log(1 + .001 * np.abs(im) ** 2), axes=0),
              extent=(fftFreq0[0], fftFreq0[-1], fftFreq1[-1], fftFreq1[0]))
    return fg

class stack_view():
    """
    Class to allow a volume to be scrubbed through with a matplotlib slider widget.
    The first axis of the volume is the slicing axis. Other keyword args are passed
    directly to imshow upon the figure creation.

    Parameters
    ----------
        stack : ndimage, 3D stack
            The stack of to show as images

    Keywords
    --------
        **kwargs :
            Passed directly to pyplot.imshow()

    """

    def __init__(self, stack, **kwargs):
        from matplotlib.widgets import Slider
        
        if stack.ndim != 3:
            raise Exception('Must be three-dimensional stack of images.')

        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        self._st = stack  # internal reference to the stack

        # Initialize the imshow axis
        self.im0 = int(self._st.shape[0] / 2)  # initial slice to show
        self.axI = self.ax.imshow(stack[self.im0, :, :], **kwargs)

        # Setup the slider
        ax_color = 'lightgoldenrodyellow'
        self.axSlider = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=ax_color)
        self.sl = Slider(self.axSlider, 'Slice', 0, self._st.shape[0] - 1, valinit=self.im0, valfmt='%1.f')
        self.sl.on_changed(self._update)

        plt.show()

    def _update(self, val):
        num = self.sl.val
        self.axI.set_data(self._st[int(round(num)), :, :])
        self.fig.canvas.draw_idle()

