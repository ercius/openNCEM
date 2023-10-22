"""
A set of visualization functions based on matplotlib.pyplot which are generally
useful for S/TEM data visualization.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import ncempy.algo
from ncempy.algo.distortion import rad_dis

def plot(dd):
    """Easy plot of data structure returned by ncempy.read(). This only works
    for 2D images currently. Calls im_calibrated with proper inputs.
    
    Parameters
    ----------
    dd : dict
        The dictionary returned by ncempy.read()
        
    Returns
    -------
    : matplotlb.pyplot.Figure
        The handle to the created Figure
    """
    
    assert dd['data'].ndim == 2
    
    return im_calibrated(dd['data'], dd['pixelSize'], units=dd['pixelUnit'])


def imsd(im, vmin=-2, vmax=2, **kwargs):
    """Show an array as an image with intensities compared to the standard deviation of the data. Other
    keyword args are passed to pyplot.imshow(). vmin and vmax are by default set to -2 and 2 respectively
    which are usually good values to set for S/TEM data.

    Parameters
    ----------
    im : np.ndarray
        The image to show.

    vmin, vmax : float, default = -2, 2
        The vmin and vmax values to pass to imshow.

    Returns
    -------
    : matplotlib.pyplot.Figure
        The handle to the created figure
    """
    fg, ax = plt.subplots(1, 1)
    im2 = im - im.mean()
    im3 = im2 / np.std(im2)
    imax = ax.imshow(im3, vmin=vmin, vmax=vmax, **kwargs)
    return fg


def im_calibrated(im, d, units=None, **kwargs):
    """ Plot an image calibrated using the pixel size d. The centers of the pixels will be the
    the center of each measurement. So, if you plot positions in real coordinates the points
    will be plotted in the center of the pixel.

    Parameters
    ---------
    im : np.ndarray
        The image to show using imshow
    d : 2-tuple
        The pixel size in both directions as a 2-tuple.
    unit : 2-tuple
        The name of the unit as a string for each pixel along each dimension. (nm or nanometer)
        
    Returns
    -------
    : pyplot.figure
        The figure containing the plot
    """
    # The default extent is calculated like this:
    ext = [-0.5, im.shape[1] - 0.5, im.shape[0] - 0.5, -0.5]
    # Calibrate the extent
    ext[0] = ext[0] * d[1]
    ext[1] = ext[1] * d[1]
    ext[2] = ext[2] * d[0]
    ext[3] = ext[3] * d[0]

    fg, ax = plt.subplots(1, 1)
    ax.imshow(im, extent=ext, **kwargs)
    
    if units:
        ax.set(xlabel=f'X ({units[1]})', ylabel=f'Y ({units[0]})')
    return fg


def imfft(im, d=1.0, ax=None, **kwargs):
    """ Show a 2D FFT as a diffractogram with log scaling applied and zero frequency
    fftshifted tp the center. A new figure is created or an axis can be specified.

    The diffracotgram is calculated from the original intensities (I) as

    .. math::
        1 + 0.001 * I ^2

    Parameters
    ----------
    im: np.ndarray
        The 2D fft of the diffraction pattern to display as a diffractogram
    d: float, optional, default = 1.0
        The real space pixel size of the image used to get the FFT
    ax: pyplot axis, optional
        An axis to plot into.
    Returns
    -------
    : matplotlib.image.Figure
        The Figure that contains the image displayed.

    Example
    -------
    This example shows how to display a 2D ndarray (image) as a
    diffractogram. The image has a real space pixel size of 0.1 nanometer.

    >> imageFFT = np.fft.fft2(im)
    >> ncempy.viz.imfft(imageFFT, d = 0.1)

    """

    fftFreq0 = np.fft.fftshift(np.fft.fftfreq(im.shape[0], d))
    fftFreq1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d))
    if ax is None:
        fg, ax = plt.subplots(1, 1)
    imax = ax.imshow(np.fft.fftshift(np.abs(im)),
                     extent=(fftFreq0[0], fftFreq0[-1], fftFreq1[-1], fftFreq1[0]), norm=colors.LogNorm(), **kwargs)
    return imax.get_figure()


def imrfft(im, d=1.0, ax=None, **kwargs):
    """Show a 2D rFFT (real FFT) as a diffractogram with log scaling applied
    and fftshift-ed along axis 0. See imfft for full details.

    Parameters
    ----------
    im : ndarray
        The 2D fft of the diffraction pattern to display as a diffractogram
    d : float, optional, default = 1.0
        The real space pixel size of the image used to get the FFT
    ax : pyplot axis, optional
        An axis to plot into.
    Returns
    -------
    : matplotlib.image.Figure
        The AxesImage that contains the image displayed.
    """

    fftFreq1 = np.fft.fftshift(np.fft.fftfreq(im.shape[1], d))
    fftFreq0 = np.fft.rfftfreq(im.shape[0], d)
    if ax is None:
        fg, ax = plt.subplots(1, 1)
    axim = ax.imshow(np.fft.fftshift(np.abs(im)), axes=0,
                     extent=(fftFreq0[0], fftFreq0[-1], fftFreq1[-1], fftFreq1[0]), norm=colors.LogNorm(), **kwargs)

    return axim.get_figure()


def im_and_fft(im, d=1.0, fft=None):
    """ Show the image and its fft side by side. Uses imfft to show the fft.

    Parameters
    ----------
    im : np.ndarray
        The image to show in both real and FFT space
    d : float
        The pixel spacing
    fft : np.ndarray, optional
        The FFT to display. If not provided then np.fft.fft2 is used.

    Returns
    -------
    : plt.figure
        The matplotlib.pyplot figure
    """
    fg, ax = plt.subplots(1, 2)
    ax[0].imshow(im, **kwargs)

    if not fft:
        fft = np.fft.fft2(im)
    imfft(fft, d=d, ax=ax[1])
    
    return fg

class stack_view:
    """
    Class to allow a volume to be scrubbed through with a matplotlib slider widget.
    The first axis of the volume is the slicing axis. Other keyword args are passed
    directly to imshow upon the figure creation.

    Parameters
    ----------
    stack : numpy.ndarray, 3D stack
        The stack of to show as images

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
    """Plot points in polar coordinate system.

    Parameters
    ----------
    points : np.ndarray
        Positions in polar coords.
    dims : tuple
        Dimension information to plot labels.
    show : bool
        Set to directly show plot in interactive mode.

    Returns
    -------
    : numpy.ndarray
        Image of the plot.

    """

    try:
        # try to convert input to np.ndarray with 2 columns (necessary if only one entry provided)
        points = np.reshape(np.array(points), (-1, 2))
        # check if enough dims available
        assert (len(dims) >= 2)
        assert (len(dims[0]) == 3)
    except:
        raise TypeError('Something wrong with the input!')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # mean value as line
    ax.axhline(np.mean(points[:, 0]), ls='--', c='k')

    # points
    ax.plot(points[:, 1], points[:, 0], 'rx')

    # labels
    ax.set_xlabel('theta /[rad]')
    ax.set_xlim((-np.pi, np.pi))
    ax.set_ylabel('r /{}'.format(dims[0][2]))

    if show:
        plt.show(block=False)

    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return plot


def plot_distpolar(points, dims, dists, ns, show=False):
    """Plot the results of distortion fitting in polar coordinates.

    Parameters
    ----------
    points : np.ndarray
        Points in polar coords.

    dims : tuple
        Dimensions, necessary to have unit information.

    dists : np.ndarray
        Results of dist fitting, length according to ns.

    ns : list
        List of used orders.

    show : bool
        Set to directly show the plot in interactive mode.

    Returns
    -------
    : np.ndarray
        Image of the plot.

    """

    try:
        # check points
        assert (isinstance(points, np.ndarray))
        assert (points.shape[1] == 2)

        # check if enough dims available
        assert (len(dims) >= 2)
        assert (len(dims[0]) == 3)

        # check orders
        assert (len(ns) >= 1)

        # check dists
        assert (dists.shape[0] == len(ns) * 2 + 1)
    except:
        raise TypeError('Something wrong with the input!')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # stuff from the single orders
    ax.axhline(dists[0], ls='--', c='k')
    xpl_ell = np.linspace(-np.pi, np.pi, 100)
    for i in range(len(ns)):
        plt.plot(xpl_ell, dists[0] * rad_dis(xpl_ell, dists[i * 2 + 1], dists[i * 2 + 2], ns[i]), 'm--')

    # points before
    ax.plot(points[:, 1], points[:, 0], 'rx')

    # sum of all distorts
    sum_dists = np.ones(xpl_ell.shape) * dists[0]
    for i in range(len(ns)):
        sum_dists *= rad_dis(xpl_ell, dists[i * 2 + 1], dists[i * 2 + 2], ns[i])
    plt.plot(xpl_ell, sum_dists, 'b-')

    # points after
    points_corr = np.copy(points)
    for i in range(len(ns)):
        points_corr[:, 0] /= rad_dis(points[:, 1], dists[i * 2 + 1], dists[i * 2 + 2], ns[i])
    plt.plot(points_corr[:, 1], points_corr[:, 0], 'gx')

    # labels
    ax.set_xlabel('theta /[rad]')
    ax.set_xlim((-np.pi, np.pi))
    ax.set_ylabel('r /{}'.format(dims[0][2]))

    if show:
        plt.show(block=False)

    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return plot


def plot_points(img, points, vminmax=(0, 1), dims=None, invert=False, show=False):
    """Plot the detected points on the input image for checking.

    Parameters
    ----------
    img : np.ndarray
        Image.
    points : np.ndarray
        Array containing the points.
    vminmax : tuple
        Tuple of two values for relative lower and upper cut off to display image.
    dims : tuple
        Tuple of dims to plot in dimensions.
    invert : bool
        Set to invert the image.
    show : bool
        Set to directly show the plot interactively.

    Returns
    -------
    : np.ndarray
        Image of the plot.

    """

    try:
        assert (isinstance(img, np.ndarray))

        assert (isinstance(points, np.ndarray))
        assert (points.shape[1] == 2)
        assert (len(points.shape) == 2)
    except:
        raise TypeError('Something wrong with the input!')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if invert:
        cmap = "Greys"
    else:
        cmap = "gray"

    if dims:
        ax.imshow(img, cmap=cmap, vmin=np.min(img) + vminmax[0] * (np.max(img) - np.min(img)),
                  vmax=np.min(img) + vminmax[1] * (np.max(img) - np.min(img)),
                  extent=(np.min(dims[0][0]), np.max(dims[0][0]), np.max(dims[1][0]), np.min(dims[1][0])))
        ax.set_xlabel('{} {}'.format(dims[0][1], dims[0][2]))
        ax.set_ylabel('{} {}'.format(dims[1][1], dims[1][2]))
        ax.set_xlim((np.min(dims[0][0]), np.max(dims[0][0])))
        ax.set_ylim((np.max(dims[1][0]), np.min(dims[1][0])))
    else:
        ax.imshow(img, cmap=cmap, vmin=np.min(img) + vminmax[0] * (np.max(img) - np.min(img)),
                  vmax=np.min(img) + vminmax[1] * (np.max(img) - np.min(img)))
        ax.set_xlim((0, img.shape[1] - 1))
        ax.set_ylim((img.shape[0] - 1, 0))

    ax.scatter(points[:, 1], points[:, 0], color='r', marker='o', facecolors='none')

    if show:
        try:
            plt.show(block=False)
        except:
            pass

    fig.canvas.draw()

    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return plot


def plot_radialprofile(r, intens, dims, show=False):
    """Plot radial profile.

    Parameters
    ----------
    r : np.ndarray
        r-axis of radial profile.
    intens : np.ndarray
        Intensity-axis of radial profile.
    dims : tuple
        Dimensions of original image to read out units.
    show : bool
        Set to directly show plot interactively.

    Returns
    -------
    : np.ndarray
        Image of the plot.
    """

    try:
        # check data
        assert (isinstance(r, np.ndarray))
        assert (isinstance(intens, np.ndarray))
        assert (np.array_equal(r.shape, intens.shape))

        # check if dims available
        assert (len(dims) >= 1)
        assert (len(dims[0]) == 3)

    except:
        raise TypeError('Something wrong with the input!')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(r, intens, 'r-')

    # labels
    ax.set_xlabel('r /{}'.format(dims[0][2]))
    ax.set_ylabel('I /[a.u.]')

    if show:
        plt.show(block=False)

    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return plot


def plot_fit(r, intens, dims, funcs, param, show=False):
    """Plot the fit results to the radial profile.

    Parameters
    ----------
    r : np.ndarray
        r-axis of radial profile.
    intens : np.ndarray
        Intensity-axis of radial profile.
    dims : tuple
        Dimensions of original image to read out units.
    funcs : tuple
        List of functions.
    param : np.ndarray
        Parameters for functions in funcs.
    show : bool
        Set to directly show plot interactively.

    Returns
    -------
    : np.ndarray
        Image of the plot.

    """

    try:
        # check data
        assert (isinstance(r, np.ndarray))
        assert (isinstance(intens, np.ndarray))
        assert (np.array_equal(r.shape, intens.shape))

        # check if dims available
        assert (len(dims) >= 1)
        assert (len(dims[0]) == 3)

        # funcs and params
        assert (len(funcs) >= 1)
        for i in range(len(funcs)):
            assert (funcs[i] in ncempy.algo.math.lkp_funcs)

        param = np.array(param)
        param = np.reshape(param, sum(map(lambda x: ncempy.algo.math.lkp_funcs[x][1], funcs)))

    except:
        raise TypeError('Something wrong with the input!')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot radial profile
    ax.plot(r, intens, 'r-')

    # plot single
    n = 0
    for i in range(len(funcs)):
        ax.plot(r, ncempy.algo.math.lkp_funcs[funcs[i]][0](r, param[n:n + ncempy.algo.math.lkp_funcs[funcs[i]][1]]),
                'g-')
        n += ncempy.algo.math.lkp_funcs[funcs[i]][1]
    # sum of functions
    sum_funcs = ncempy.algo.math.sum_functions(r, funcs, param)
    ax.plot(r, sum_funcs, 'b-')

    # labels
    ax.set_xlabel('r /{}'.format(dims[0][2]))
    ax.set_ylabel('I /[a.u.]')

    if show:
        plt.show(block=False)

    # render to array
    fig.canvas.draw()
    plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return plot
