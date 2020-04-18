"""
A set of visualization functions based on matplotlib.pyplot which are generally
useful for S/TEM data visualization.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt

from ncempy.algo.distortion import rad_dis


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


def plot_ringpolar(points, dims, show=False):
    '''Plot points in polar coordinate system.

    Parameters:
        points (np.ndarray):    Positions in polar coords.
        dims (tuple):    Dimension information to plot labels.
        show (bool):    Set to directly show plot in interactive mode.

    Returns:
        (np.ndarray):    Image of the plot.

    '''

    try:
        # try to convert input to np.ndarray with 2 columns (necessary if only one entry provided)
        points = np.reshape(np.array(points), (-1,2))
        # check if enough dims available
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


def plot_distpolar(points, dims, dists, ns, show=False):
    '''Plot the results of distortion fitting in polar coordinates.

    Parameters:
        points (np.ndarray):    Points in polar coords.
        dims (tuple):    Dimensions, necessary to have unit information.
        dists (np.ndarray):    Results of dist fitting, length according to ns.
        ns (list):    List of used orders.
        show (bool):    Set to directly show the plot in interactive mode.

    Returns:
        (np.ndarray):    Image of the plot.

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