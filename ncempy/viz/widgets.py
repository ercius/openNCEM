"""ipywidgets-based interactive visualization tools for use in Jupyter notebooks.

Requires the 'jupyter' extra::

    pip install ncempy[jupyter]

Functions
---------
SI_interactive      : Browse a spectrum image by spatial position.
image_and_plot      : Scroll through an image series with a linked line plot.
fourDstem_interactive : Browse 4D-STEM diffraction patterns by probe position.
stack_browse        : Scroll through a 3D image stack.
stemtomo_browse     : Browse STEMTomo7 data by tilt angle and image index.
bandpass_interactive : Interactively apply a bandpass filter to an image.
volume_slicer       : View orthogonal slices through a 3D volume.
"""

try:
    import ipywidgets as widgets
    from ipywidgets import interactive
except ImportError:
    raise ImportError(
        "ipywidgets is required for ncempy.viz.widgets. "
        "Install it with: pip install ncempy[jupyter]"
    )

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def SI_interactive(data, box_size):
    """Browse a spectrum image by moving a spatial position selector.

    Parameters
    ----------
    data : ndarray, shape (energy, Y, X)
        The spectrum image dataset.
    box_size : int
        Half-width of the summation box around the selected position. Must be >= 1.

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    if box_size < 1:
        raise ValueError("box_size must be >= 1")

    b_2 = 1 if box_size == 1 else box_size // 2

    im = data.sum(axis=0).copy()
    spec = data[:, 0:box_size, 0:box_size].sum(axis=(1, 2))

    fg, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 4))
    p0 = ax2.plot(spec)
    ax2.set_yscale('log')
    imax1 = ax1.imshow(im)
    p1 = ax1.scatter(b_2, b_2, c='r')
    fg.tight_layout()

    def axUpdate(ii, jj):
        new_spec = data[:, ii - b_2:ii + b_2, jj - b_2:jj + b_2].sum(axis=(1, 2))
        p0[0].set_ydata(new_spec)
        ax2.set_ylim(new_spec.min(), new_spec.max())
        p1.set_offsets((jj, ii))

    w1 = widgets.IntSlider(
        value=b_2, min=b_2, max=data.shape[1] - b_2,
        description='Y:', continuous_update=True,
    )
    w2 = widgets.IntSlider(
        value=b_2, min=b_2, max=data.shape[2] - b_2,
        description='X:', continuous_update=True,
    )
    return interactive(axUpdate, ii=w1, jj=w2)


def image_and_plot(images, plots):
    """Scroll through an image series with a linked line plot.

    Parameters
    ----------
    images : ndarray, shape (N, Y, X)
        Stack of images to display.
    plots : ndarray, shape (N, M)
        Corresponding 1D data to plot for each image (mean-subtracted on display).

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    fg, ax = plt.subplots(1, 2)
    imax = ax[0].imshow(images[0])
    p1 = ax[1].plot(plots[0, :])

    def ax_update(ii):
        imax.set_array(images[ii])
        y_data = plots[ii, :] - plots[ii, :].mean()
        p1[0].set_ydata(y_data)
        ax[1].set_ylim(y_data.min(), y_data.max())

    w1 = widgets.IntSlider(
        value=0, min=0, max=len(images) - 1,
        description='index:', continuous_update=True, readout=False,
    )
    return interactive(ax_update, ii=w1)


def fourDstem_interactive(data, image):
    """Browse 4D-STEM diffraction patterns by probe position.

    Parameters
    ----------
    data : ndarray, shape (Y, X, ky, kx)
        The 4D-STEM dataset.
    image : ndarray, shape (Y, X)
        Real-space image (e.g. bright-field) used to select the probe position.

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    fg, (ax1, ax2) = plt.subplots(1, 2)
    imax1 = ax1.imshow(image)
    p1 = ax1.plot(image.shape[1] // 2, image.shape[0] // 2, 'or')[0]
    imax2 = ax2.imshow(data[data.shape[0] // 2, data.shape[1] // 2, :, :], norm=LogNorm())
    ax1.set(title='Real-space image', xlabel='X', ylabel='Y')
    ax2.set(title='Diffraction pattern (log scale)')
    fg.tight_layout()

    def axUpdate(j, i):
        p1.set_xdata((i,))
        p1.set_ydata((j,))
        dp = data[j, i, :, :]
        imax2.set_array(dp)
        imax2.set_clim(dp.min(), dp.max())

    w_j = widgets.IntSlider(
        value=image.shape[0] // 2, min=0, max=image.shape[0] - 1,
        description='Y:', continuous_update=True,
    )
    w_i = widgets.IntSlider(
        value=image.shape[1] // 2, min=0, max=image.shape[1] - 1,
        description='X:', continuous_update=True,
    )
    return interactive(axUpdate, j=w_j, i=w_i)


def stack_browse(data, **kwargs):
    """Scroll through a 3D image stack slice by slice.

    Parameters
    ----------
    data : ndarray, shape (N, Y, X)
        The image stack. The first axis is the slicing axis.
    **kwargs
        Passed to pyplot.imshow().

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    fg, ax = plt.subplots(1, 1, figsize=(4, 4))
    imax = ax.imshow(data[0, :, :], **kwargs)
    fg.tight_layout()

    def axUpdate(ii):
        tt = data[ii, :, :]
        imax.set_data(tt)
        imax.set_clim(tt.min(), tt.max())

    w1 = widgets.IntSlider(
        value=0, min=0, max=data.shape[0] - 1,
        description='slice:', continuous_update=True,
    )
    return interactive(axUpdate, ii=w1)


def stemtomo_browse(data, **kwargs):
    """Browse STEMTomo7 data by tilt angle and image index.

    Parameters
    ----------
    data : ndarray, shape (tilts, images, Y, X)
        Data output from the STEMTomo7 acquisition program.
    **kwargs
        Passed to pyplot.imshow().

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    fg, ax = plt.subplots(1, 1, figsize=(4, 4))
    imax = ax.imshow(data[0, 0, :, :], **kwargs)
    fg.tight_layout()

    def axUpdate(ii, jj):
        tt = data[ii, jj, :, :]
        imax.set_data(tt)
        imax.set_clim(tt.min(), tt.max())

    w1 = widgets.IntSlider(
        value=0, min=0, max=data.shape[0] - 1,
        description='tilt:', continuous_update=True,
    )
    w2 = widgets.IntSlider(
        value=0, min=0, max=data.shape[1] - 1,
        description='image #:', continuous_update=True,
    )
    return interactive(axUpdate, ii=w1, jj=w2)


def bandpass_interactive(im):
    """Interactively apply a bandpass filter to a 2D image.

    Shows the filtered image and its FFT magnitude side by side. The inner and
    outer frequency cutoffs and the Gaussian smoothing width are adjustable.

    Parameters
    ----------
    im : ndarray, shape (Y, X)
        The input image to filter.

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    def _bandpass_filter(im, inner, outer, sigma):
        num0, num1 = im.shape
        nf0 = np.fft.fftshift(np.fft.fftfreq(num0))
        nf1 = np.fft.fftshift(np.fft.fftfreq(num1))
        XX, YY = np.meshgrid(nf0, nf1, indexing='ij')
        mask = (np.sqrt(XX**2 + YY**2) < outer) & (np.sqrt(XX**2 + YY**2) > inner)

        linG0 = np.linspace(-(num0 / 2), num0 / 2 - 1, num0)
        linG1 = np.linspace(-(num1 / 2), num1 / 2 - 1, num1)
        gXX, gYY = np.meshgrid(linG0, linG1, indexing='ij')
        gg = np.exp(-(gXX**2 + gYY**2) / (2 * sigma**2)) / (2 * np.pi * sigma**2)

        mask_smooth = np.fft.irfft2(np.fft.rfft2(gg) * np.fft.rfft2(mask))
        return np.real(np.fft.ifft2(np.fft.fft2(im) * mask_smooth))

    fg, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    imax1 = ax1.imshow(im)
    imax2 = ax2.imshow(np.abs(np.fft.fftshift(np.fft.rfft2(im), axes=0)))
    ax1.set(title='Filtered image')
    ax2.set(title='FFT magnitude')
    fg.tight_layout()

    def axUpdate(io, sigma):
        tt = _bandpass_filter(im, io[0], io[1], sigma)
        tt_fft = np.abs(np.fft.fftshift(np.fft.rfft2(tt), axes=0))
        imax1.set_data(tt)
        imax1.set_clim(tt.min(), tt.max())
        imax2.set_data(tt_fft)
        imax2.set_clim(tt_fft.min(), tt_fft.max())

    w_bandpass = widgets.FloatRangeSlider(
        value=[0.1, 0.4], min=0, max=0.5, step=0.005,
        description='bandpass:', continuous_update=True, readout_format='.3f',
    )
    w_sigma = widgets.FloatSlider(
        value=2, min=0.1, max=5, step=0.1,
        description='sigma:', continuous_update=True, readout_format='.1f',
    )
    return interactive(axUpdate, io=w_bandpass, sigma=w_sigma)


def volume_slicer(data):
    """View orthogonal slices through a 3D volume with linked crosshair markers.

    Displays XY, XZ, and YZ slices with a red crosshair showing the current
    slice position, plus a summed projection for reference.

    Parameters
    ----------
    data : ndarray, shape (Z, Y, X)
        The 3D volume to explore.

    Returns
    -------
    : ipywidgets.interactive
        The interactive widget. Call display() on the result to show it.
    """
    s = [ii // 2 for ii in data.shape]
    vmin, vmax = data.min(), data.max()

    fg, ax = plt.subplots(2, 2)
    imax00 = ax[0, 0].imshow(data[s[0], :, :], vmin=vmin, vmax=vmax)
    p00 = ax[0, 0].plot(s[2], s[1], '+r')
    ax[0, 0].set(title='XY slice')

    imax10 = ax[1, 0].imshow(data[:, s[1], :], vmin=vmin, vmax=vmax)
    p10 = ax[1, 0].plot(s[2], s[0], '+r')
    ax[1, 0].set(title='XZ slice')

    imax11 = ax[1, 1].imshow(data[:, :, s[2]], vmin=vmin, vmax=vmax)
    p11 = ax[1, 1].plot(s[1], s[0], '+r')
    ax[1, 1].set(title='YZ slice')

    ax[0, 1].imshow(data.sum(axis=1), vmin=vmin * data.shape[1], vmax=vmax * data.shape[1])
    ax[0, 1].set(title='Y projection')
    fg.tight_layout()

    def axUpdate(i, j, k):
        imax00.set_data(data[i, :, :])
        p00[0].set_xdata((k,))
        p00[0].set_ydata((j,))
        imax10.set_data(data[:, j, :])
        p10[0].set_xdata((k,))
        p10[0].set_ydata((i,))
        imax11.set_data(data[:, :, k])
        p11[0].set_xdata((j,))
        p11[0].set_ydata((i,))

    return interactive(
        axUpdate,
        i=(0, data.shape[0] - 1),
        j=(0, data.shape[1] - 1),
        k=(0, data.shape[2] - 1),
    )
