"""
Module to find local maxima in 2D and 3D data such as images and density maps.

Many functions are based on work by Colin Ophus, cophus@lbl.gov and his excellent
RealSpaceLattice01.m code.

author: Peter Ercius, percius@lbl.gov
"""

try:
    from . import gaussND
except:
    try:
        print('peakFind: Temporary import!')
        import gaussND
    except ModuleNotFoundError:
        raise ModuleNotFoundError
    

import numpy as np
import scipy.optimize as opt
from scipy import ndimage

default_indexing: str = 'ij'


def doubleRoll(image, vec):
    """ Roll an 2D ndarray (e.g. an image) along each dimension. 
    
        Parameters
        ----------
            image : np.ndarray
                The 2D array to roll.
            
            vec : iterable
                The number of time to apply the roll to the rows and columns where
                the tuple has the structure (row, col). Values in tuple must be int.
                
        Returns
        -------
            : np.ndarray
                The rolled image.
    """
    assert len(vec) == 2
    return np.roll(np.roll(image, vec[0], axis=0), vec[1], axis=1)


def tripleRoll(vol, vec):
    """ Roll a 3D ndarray (e.g. an volume) along each dimension. 
    
        Parameters
        ----------
            vol : np.ndarray
                The 3D array to roll
            
            vec : iterable
                The number of time to apply the roll the tuple has the structure (axis0, axis1, axis2).
                Values in tuple must be int.
                
        Returns
        -------
            : ndarray
                The rolled volume.
    """
    assert len(vec) == 3
    return np.roll(np.roll(np.roll(vol, vec[0], axis=0), vec[1], axis=1), vec[2], axis=2)


def peakFind2D(image, threshold):
    """ Find peaks in a 2D image based on the roll technique. A threshold is first applied to remove
    low intensity noisy peaks.
    
    Parameters
    ----------
        image : ndarray
            An 2D ndarray of image intensities with peaks to find.
        
        threshold : float
            Threshold [0,1) the image to remove noisy peaks based on max intensity in the image
        
    Returns
    -------
        positions : array
            A Nx2 array with the index locations of the n number of peaks
    """
    if (threshold < 0) or (threshold >= 1):
        print('Error: Threshold must be 0 <= threshold < 1')
        return 0

    # Could generate the set of coordinates using
    # list(set(permutations((0,1,1,-1,-1),2)))
    pLarge = (image > doubleRoll(image, [-1, -1])) & \
             (image > doubleRoll(image, [0, -1])) & \
             (image > doubleRoll(image, [1, -1])) & \
             (image > doubleRoll(image, [-1, 0])) & \
             (image > doubleRoll(image, [1, 0])) & \
             (image > doubleRoll(image, [-1, 1])) & \
             (image > doubleRoll(image, [0, 1])) & \
             (image > doubleRoll(image, [1, 1])) & \
             (image > threshold * np.max(image))

    positions = np.array(np.where(pLarge * image)).T.copy()  # return array as [num,posXY]
    return positions.astype(int)


def peakFind3D(vol, threshold):
    """
    Find peaks in a 3D volume based on the roll technique. A threshold is first applied to remove
    low intensity noisy peaks.
    
    Parameters
    ------
        vol : ndarray
            An 3D ndarray of image intensities with peaks to find.

        threshold : float
            Threshold [0,1) the image to remove noisy peaks based on max intensity in the image
        
    Returns
    -------
        positions : ndarray
            A Nx3 array with the index locations of the n number of peaks
    """
    if (threshold < 0) or (threshold >= 1):
        print('Error: Threshold must be 0 <= threshold < 1')
        return 0

    # Could generate the set of coordinates programmatically
    # using list(set(permutations((0,0,1,1,1,-1,-1,-1),3))
    pLarge = ((vol > tripleRoll(vol, [-1, -1, -1]))
              & (vol > tripleRoll(vol, [0, -1, -1]))
              & (vol > tripleRoll(vol, [1, -1, -1]))
              & (vol > tripleRoll(vol, [-1, 0, -1]))
              & (vol > tripleRoll(vol, [1, 0, -1]))
              & (vol > tripleRoll(vol, [-1, 1, -1]))
              & (vol > tripleRoll(vol, [0, 1, -1]))
              & (vol > tripleRoll(vol, [1, 1, -1]))
              & (vol > tripleRoll(vol, [0, 0, -1]))
              & (vol > tripleRoll(vol, [-1, -1, 0]))
              & (vol > tripleRoll(vol, [0, -1, 0]))
              & (vol > tripleRoll(vol, [1, -1, 0]))
              & (vol > tripleRoll(vol, [-1, 0, 0]))
              & (vol > tripleRoll(vol, [1, 0, 0]))
              & (vol > tripleRoll(vol, [-1, 1, 0]))
              & (vol > tripleRoll(vol, [0, 1, 0]))
              & (vol > tripleRoll(vol, [1, 1, 0]))
              & (vol > tripleRoll(vol, [-1, -1, 1]))
              & (vol > tripleRoll(vol, [0, -1, 1]))
              & (vol > tripleRoll(vol, [1, -1, 1]))
              & (vol > tripleRoll(vol, [-1, 0, 1]))
              & (vol > tripleRoll(vol, [1, 0, 1]))
              & (vol > tripleRoll(vol, [-1, 1, 1]))
              & (vol > tripleRoll(vol, [0, 1, 1]))
              & (vol > tripleRoll(vol, [1, 1, 1]))
              & (vol > tripleRoll(vol, [0, 0, 1]))
              & (vol > threshold * np.max(vol)))

    positions = np.array(np.where(pLarge * vol)).T.copy()  # return array as [num,posXYZ]
    return positions.astype(int)


def enforceMinDist(positions, intensities, minDistance):
    """ Remove peaks that violate a minimum distance requirement. If two peaks are
    closer together than the minDistance then the one with the highest intensity is chosen.
    
    Works for 2D and 3D data sets.
    
    Parameters
    ----------
        positions : np.ndarray
            An array of shape Nx3 or Nx2 containing the positions of the peaks to test 
        intensities : np.ndarray
            An array of shape Nx1 with the peak intensities.
        minDistance : int
            The minimum distance in pixels
        
    Returns
    -------
        validPeaks : np.ndarray
            An ndarray of Mx3 or Mx2 (M <= N) with invalid peaks removed.
    """

    validPeaks = np.zeros(positions.shape, dtype=np.float32)  # Will hold intensity,X,Y,Z

    alreadyTested = np.zeros(
        positions.shape[0])  # all peaks are tested unless they were already compared to another peak

    for ii, p0 in enumerate(positions):  # ii is the iteration number and p0 is the element
        if not alreadyTested[ii]:
            pointsDist = np.sqrt(np.sum((positions - p0) ** 2, axis=1))  # find the distance to every other point
            closePeaks = pointsDist[(pointsDist > 0) & (pointsDist < minDistance)]

            if closePeaks.any():  # test for any values < minDistance
                closePeaksIndices = np.where(pointsDist < minDistance)[0]  # find all close peaks including self
                peakIntensities = intensities[closePeaksIndices]
                highestPeakIndex = closePeaksIndices[np.argmax(peakIntensities)]  # find the index of the highest peak
                alreadyTested[closePeaksIndices] = 1  # Don't test peaks in the future if already tested
                validPeaks[closePeaksIndices, :] = np.nan  # set all values to nan; cant be integer data set
                validPeaks[ii, :] = positions[highestPeakIndex, :]  # set highest intensity to correct values
            else:
                validPeaks[ii, :] = positions[ii, :]  # save the current peak
                # alreadyTested[ii] = 1 #Dont test peaks in the future if already tested

    validPeaks = validPeaks[~np.isnan(validPeaks[:, 0]), :]
    # print("Removed %d peaks" % (positions.shape[0] - validPeaks.shape[0]))
    return validPeaks.astype(int)


def peaksToVolume(peakList, volShape, gaussSigma, gaussSize, indexing='ij'):
    """ Place 3D Gaussian at a set of peak positions.
    
    Parameters
    ----------
        peakList : np.ndarray
            Array of peak positions of size (numPeaks, 3).
        
        volShape : tuple
            3 element tuple of the shape of the volume containing the peaks
        
        gaussSigma : tuple
            3 element tuple with that sigma parameters passed to the gaussND.gauss3D() function (sigX,sigY,sigZ)
        
        gaussSize : tuple
            3 element tuple describing the 3D size of the simulated gaussian box. Must be odd.

        indexing : str
            String to pass to np.meshgrid to indicate which indexing is used. Either 'ij' or 'xy' are possible. Default
            is 'ij' for this module.
    Returns
    -------
        : np.ndarray
            A volume containing gaussian peaks at each position indicated in peakList parameter.
    
    """

    sim_volume = np.zeros(volShape, dtype=np.float32)

    if np.any(np.array(gaussSize) % 2 == 0):
        raise ValueError('gaussSize must be odd.')

    sz = [int((g - 1) // 2) for g in gaussSize]
    temp0 = np.arange(-sz[0], sz[0]+1)
    temp1 = np.arange(-sz[1], sz[1]+1)
    temp2 = np.arange(-sz[2], sz[2]+1)
    Y3D, X3D, Z3D = np.meshgrid(temp0, temp1, temp2, indexing=indexing)

    for pos in peakList:
        pos_loc = (pos // 1).astype(int)

        pos_remain = pos % 1
        gg3 = gaussND.gauss3D(Y3D, X3D, Z3D,
                              pos_remain[0], pos_remain[1], pos_remain[2],
                              gaussSigma[0], gaussSigma[1], gaussSigma[2])

        sim_volume[pos_loc[0] - sz[0]:pos_loc[0] + sz[0] + 1,
                   pos_loc[1] - sz[1]:pos_loc[1] + sz[1] + 1,
                   pos_loc[2] - sz[2]:pos_loc[2] + sz[2] + 1] += gg3
    return sim_volume


def peaksToImage(peakList, imShape, gaussSigma, gaussSize, indexing='ij'):
    """ Place 2D Gaussian at a set of peak positions.
    
    Parameters
    ----------
        peakList : np.ndarray
            Array of peak positions of size (numPeaks, 2).
        
        imShape : tuple
            2 element tuple of the shape of the image containing the peaks
        
        gaussSigma : tuple
            2 element tuple with that sigma parameters passed to the gaussND.gauss2D() function (sigX,sigY)
        
        gaussSize : tuple
            2 element tuple describing the 2D size of the simulated gaussian box. Must be odd.

        indexing : str, default = 'ij'
            'ij' or 'xy' indexing to pass to np.meshgrid.
    Returns
    -------
        : np.ndarray
            A image containing Gaussian peaks at each position indicated in peakList parameter.
    
    """

    sim_image = np.zeros(imShape, dtype=np.float32)

    if np.any(np.array(gaussSize) % 2 == 0):
        raise ValueError('gaussSize must be odd.')

    sz = [int((g - 1) // 2) for g in gaussSize]
    temp0 = np.arange(-sz[0], sz[0] + 1)
    temp1 = np.arange(-sz[1], sz[1] + 1)
    Y3D, X3D = np.meshgrid(temp1, temp0, indexing=indexing)
    for pos in peakList:
        pos_loc = (pos // 1).astype(int)
        pos_remain = pos % 1
        gg2 = gaussND.gauss2D(Y3D, X3D,
                              pos_remain[0], pos_remain[1],
                              gaussSigma[0], gaussSigma[1])
        sim_image[pos_loc[0] - sz[0]:pos_loc[0] + sz[0] + 1,
                  pos_loc[1] - sz[1]:pos_loc[1] + sz[1] + 1] += gg2
    return sim_image


def lattice2D_norm(u, v, a, b, origin, num_points):
    """
    Returns a set of points in a lattice according to the u, v unit vectors (vectors are normalized internally)
    and lengths a,b centered at origin. The lattice has num_points along each u,v
    vector.
    
    Parameters
    ----------
        u, v : 2-tuple
            2 element tuples defining the lattice directions as vectors.

        a, b : float
            values to multiply each vector by (if u,v are not normalized then set these to 1)

        origin : 2-tuple
            The origin in the format (x0, y0)

        num_points : 2-tuple
            2 element tuple of the number of repeats of the lattice along each direction

    Returns
    ------
        : ndarray
            The set of points in the lattice. Size (num_points[0] * num_points[1], 2)
    """

    assert np.linalg.norm(u) == 1
    assert np.linalg.norm(v) == 1

    totalNumPoints = np.prod(num_points)

    Y2D, X2D = np.meshgrid(range(num_points[1]), range(num_points[0]), indexing=default_indexing)

    X = X2D.reshape(totalNumPoints)
    Y = Y2D.reshape(totalNumPoints)

    xy = np.zeros((np.prod(num_points), 2))
    xy[:, 0] = origin[0] + a * X * u[0] + b * Y * v[0]
    xy[:, 1] = origin[1] + a * X * u[1] + b * Y * v[1]

    return xy


def lattice2D(u, v, a, b, origin, num_points):
    """
    A modified version of peakFind.lattice2D_norm which can use non normalized u,v,vectors
    
    Parameters
    ----------
        u, v : tuple
            2 element vectors defining the lattice directions

        a, b : float
            values to multiply each vector by (if u,v,w are not normalized then set these to 1)

        origin : tuple
            The origin in the format (x0,y0)

        num_points : tuple
            2 element tuple of the number of repeats of the lattice along each direction

    Returns
    ------
        : ndarray
            The set of points in the lattice. Size (num_points[0] * num_points[1], 2)
    """
    totalNumPoints = np.prod(num_points)

    Y2D, X2D = np.meshgrid(range(num_points[1]), range(num_points[0]), indexing=default_indexing)

    X = X2D.reshape(totalNumPoints)
    Y = Y2D.reshape(totalNumPoints)

    xy = np.zeros((np.prod(num_points), 2))
    xy[:, 0] = origin[0] + a * X * u[0] + b * Y * v[0]
    xy[:, 1] = origin[1] + a * X * u[1] + b * Y * v[1]

    return xy


def lattice3D(u, v, w, a, b, c, origin, num_points):
    """
    Returns a set of points in a lattice according to the u, v, w vectors
    and lengths a,b centered at origin. The lattice has num_points along each u,v,w
    vector.

    Parameters
    ----------
        u, v, w : tuple or np.ndarray
            3 element vectors defining the lattice directions
        a, b, c : float
            Values to multiply each vector by (if u,v,w are not normalized then set these to 1)
        origin : tuple
            The origin in the format (x0,y0,z0)
        num_points : tuple
            3 element tuple of the number of repeats of the lattice along each direction

    Returns
    -------
        xyz : ndarray
            A (N,3) shaped set of coordinates
    """

    totalNumPoints = np.prod(num_points)

    Y3D, X3D, Z3D = np.meshgrid(range(num_points[1]), range(num_points[0]), range(num_points[2]),
                                indexing=default_indexing)

    X = X3D.reshape(totalNumPoints)
    Y = Y3D.reshape(totalNumPoints)
    Z = Z3D.reshape(totalNumPoints)

    xyz = np.zeros((np.prod(num_points), 3))
    xyz[:, 0] = origin[0] + a * X * u[0] + b * Y * v[0] + c * Z * w[0]
    xyz[:, 1] = origin[1] + a * X * u[1] + b * Y * v[1] + c * Z * w[1]
    xyz[:, 2] = origin[2] + a * X * u[2] + b * Y * v[2] + c * Z * w[2]

    return xyz


def applyLatticeLimit(lattice, bounds):
    """Remove lattice points outside the data bounds. For 2D and 3D data.
    
    Parameters
    ---------
        lattice : ndarray; (N, 2) or (N, 3)
            From lattice2D

        bounds : tuple, 
            Minimum and maximum for axes 0 and 1
            (min0, max0, min1, max1) or axes 0, 1 and 2
            (min0, max0, min1, max1, min2, max2)

    Returns
    -------
        : ndarray; (M, 2) or (M, 3)
            Same as lattice input except only containing
            points within the bounds specified. M <= N
    """
    if len(bounds) == 4:
        goodUVs = (((lattice[:, 0] > bounds[0]) & (lattice[:, 0] < bounds[1])) &
                   ((lattice[:, 1] > bounds[2]) & (lattice[:, 1] < bounds[3])))
    elif len(bounds) == 6:
        goodUVs = (((lattice[:, 0] > bounds[0]) & (lattice[:, 0] < bounds[1])) &
                   ((lattice[:, 1] > bounds[2]) & (lattice[:, 1] < bounds[3])) &
                   ((lattice[:, 2] > bounds[4]) & (lattice[:, 2] < bounds[5])))
    else:
        print('Bounds needs be be either 4 or 6 value tuple.')
        return None
    return lattice[goodUVs, :]


# noinspection PyArgumentList
def peakPlot3D(X, Y, Z, mkr, myAxes3D):
    """
    Plot a set of peaks in a 3D plot using matplotlib. See the example below for how to set up a figure as input
    to this function using matplotlib Axes3D.

    todo: move this to viz

    Parameters
    -----------
        X : np.ndarray
            The x positions as a 1D array.

        Y : np.ndarray
            The y positions as a 1D array.

        Z : np.ndarray
            The Z positions as a 1D array.

        mkr : string
            String indicating the marks type and color. Passed directly to pyplot.plt
            command (ex. 'go' for green circles)

        myAxes3D : mpl_toolkits.mplot3d.Axes3d
            An Axes3D object. (see example below)
        
    Returns
    -------
        None
    
    Example
    -------
    This function requires the input of a Axes3D to plot a set of peaks. See below how to import
    and set up a figure for use with this function.
    >> import matplotlib.pyplot as plt
    >> import mpl_toolkits.mplot3d
    >> from ncempy.algo import peak_find
    >> fg1 = plt.figure()
    >> ax1 = mpl_toolkits.mplot3d.Axes3D(fg1)
    >> peak_find.peakPlot3D(peakList[:,2], peakList[:,1], peakList[:,0], 'go', ax1)
    """

    import matplotlib.pyplot as plt

    plt.plot(X, Y, Z, mkr)

    max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max() / 2.0

    mid_x = (X.max() + X.min()) * 0.5
    mid_y = (Y.max() + Y.min()) * 0.5
    mid_z = (Z.max() + Z.min()) * 0.5

    myAxes3D.set_xlim(mid_x - max_range, mid_x + max_range)
    myAxes3D.set_ylim(mid_y - max_range, mid_y + max_range)
    myAxes3D.set_zlim(mid_z - max_range, mid_z + max_range)

    myAxes3D.set_xlabel('X')
    myAxes3D.set_ylabel('Y')
    myAxes3D.set_zlabel('Z')


def peak3View(fg, vol, peak_positions):
    """
    Plot 3 orthogonal slices and the corresponding peaks in each slice.
    
    Not fully test. Unsure whether the slices and peaks are exactly the same.

    todo: move this to viz

    Parameters
    -----------
        fg : matplotlib figure
            The figure to use. This allows the user to reuse figures

        vol : ndarray (3D)
            The volume to slice along each orthogonal direction. Set up for [X,Z,Y] orientation

        peak_positions : ndarray, list, tuple (1D)
            This selects the point shown in all three slices
        
    Returns
    -------
        None

    """
    ((ax1, ax2), (ax3, ax4)) = fg.subplots(2, 2)

    slices = [int(round(ii)) for ii in peak_positions]
    ax1.imshow(vol[slices[0], :, :])
    ax1.plot(peak_positions[2], peak_positions[1], 'xk')
    ax1.set(title='X slice {}'.format(slices[0]))

    ax3.imshow(vol[:, slices[1], :])
    ax3.plot(peak_positions[2], peak_positions[0], 'xk')
    ax3.set(title='Y slice {}'.format(slices[1]))

    ax4.imshow(vol[:, :, slices[2]])
    ax4.plot(peak_positions[1], peak_positions[0], 'xk')
    ax4.set(title='Z slice {}'.format(slices[2]))

    ax2.axis('off')
    fg.tight_layout()


def remove_xrays(imageOriginal, threshold, size_median_filter=(3, 3)):
    """
    Find x-ray pixels in TEM images and replace the xray intensity 
    with a local median value.
    
    Parameters
    ----------
        imageOriginal : ndarray
            The image as a 2D ndarray.
        threshold : float
            Value above which to consider a pixel intensity an x-ray.
        size_median_filter : list or tuple
            The footprint of the median filter over which to search for a median value.
            
    Returns
    -------
        : ndarray
            The image with x-ray values replaced with median local value.
    
    """
    # Find x-ray pixels
    image_med_filter = ndimage.filters.median_filter(imageOriginal, size_median_filter)
    sub_xrays = np.where(imageOriginal > (image_med_filter + threshold))

    # Replace x-ray pixels with median filtered values.
    image_filtered = imageOriginal.copy()
    image_filtered[sub_xrays[0], sub_xrays[1]] = image_med_filter[sub_xrays[0], sub_xrays[1]]
    return image_filtered, sub_xrays


def writeXYZ(filename, XYZ, element, comment):
    """
    Write out a set of XYZ coordinates that can be read by various crystal viewing software such as Vesta.

    todo: Move this to io

    Parameters
    ----------
        filename : str
            The name of the file to write to.

        XYZ : np.ndarray
            Atom coordinates in a 2D ndarray with axes [atom number, position]. Should be shape [N, 3]

        element : tuple
            Element symbol of each coordinate as a string (e.g. 'C'). Should be length N.

        comment : str
            A comment to add to the file.

    Example
    -------
        # Write out the peak positions as Carbon atoms:
        >> writeXYZ('name.xyz', positions, ['C',] * len(positions), 'writeXYZ example carbon atoms')
    """
    numAtoms = XYZ.shape[0]
    with open(filename, 'w') as f1:
        f1.write('{}\n'.format(numAtoms))
        f1.write('{}\n'.format(comment))
        for ii in range(numAtoms):
            f1.write('{} {} {} {} \n'.format(element[ii], XYZ[ii, 0], XYZ[ii, 1], XYZ[ii, 2]))


def refineLattice2D(or0, u0, v0, pos, fraction=(1, 1), max_iter=30,
                    refine_locally=True, verbose=False, num_unit_cells=(1, 1)):
    """ Refine lattice based on measurements and initial
    guess at lattice vectors. This code is designed to work only with
    square lattices. Hexagonal and more complex lattice might not work.

    Parameters
    ----------
        or0 : tuple
            Initial guess of origin. 2-tuple
        u0 : tuple
            Initial u vector. Use num_unit_cells to average over several lengths of the lattice vector in this
            direction. 2-tuple
        v0 : tuple
            Initial v vector. Use num_unit_cells to average over several lengths of the lattice vector in this
            direction. 2-tuple
        pos : np.ndarray
            The set of positions to fit the lattice to. (number, position) (M, 2)
        refine_locally : bool
            Refine locally near the origin before using all positions
            Locally is considered 2 time the larger vector of u0 or v0.
        fraction : tuple
            Site fraction to take into account peak positions inside the unit cell. FCC imaged along [100] for\
            example would need site_fraction = (2, 2).
        max_iter : int
            The maximum number of iterations to run to refine. This usually converges in a few iterations. Default is
            30.
        num_unit_cells : tuple
            The number of unit cells the initial guess of u0 and v0 extend over. 2-tuple (num_u0, num_v0) with default
        verbose : bool
            Print out helpful information.
            
    Returns
    -------
        : tuple
            Tuple of of 4 np.ndarray (origin, u, v, ab) where ab is the
            site positions in fractions of u and v.
    
    """

    origin = or0
    u = u0
    v = v0

    xyRes = (0, 0)
    ab = (0, 0)

    # Divide by number plotted
    u = [ii / num_unit_cells[0] for ii in u]
    v = [ii / num_unit_cells[1] for ii in v]

    for ii in range(max_iter):
        if refine_locally is True:
            for jj in range(1, 30):
                # Refine locally first stepping out
                rMax = jj * np.mean((np.linalg.norm(u), np.linalg.norm(v)))  # radius to use first

                # compute a, b, c values by finding points p0 in basis {u,v,w} from {i,j,k}
                uv = np.array([u, v]).T
                pR = np.sqrt((pos[:, 0] - origin[0]) ** 2 + (pos[:, 1] - origin[1]) ** 2)
                p1 = pos[pR < rMax, :]
                if p1.shape[0] < 10:
                    # Test that there are enough atoms to do a fit. Skip if less than 10
                    continue
                (ab, abRes, abX, abY) = np.linalg.lstsq(uv, (p1 - origin).T, rcond=None)

                # Refine lattice
                A = np.ones([ab.shape[1], 3])  # origin values should be 1
                A[:, 1] = np.round(fraction[0] * ab[0, :]) / fraction[
                    0]  # fill the rest of the values with the abc values
                A[:, 2] = np.round(fraction[1] * ab[1, :]) / fraction[1]

                (xy_beta, xyRes, rXY, sXY) = np.linalg.lstsq(A, p1, rcond=None)

                origin = xy_beta[0, :]
                u = xy_beta[1, :]
                v = xy_beta[2, :]

        uv = np.array([u, v]).T
        (ab, abRes, abX, abY) = np.linalg.lstsq(uv, (pos - origin).T, rcond=None)

        # Refine lattice
        A = np.ones([ab.shape[1], 3])
        A[:, 1] = np.round(fraction[0] * ab[0, :]) / fraction[0]  # fill the rest of the values with the abc values
        A[:, 2] = np.round(fraction[1] * ab[1, :]) / fraction[1]
        (xy_beta, xyRes, rXY, sXY) = np.linalg.lstsq(A, pos, rcond=None)

        origin = xy_beta[0, :]
        u = xy_beta[1, :]
        v = xy_beta[2, :]
        if verbose:
            print('[u,v] = [[{}],[{}]]'.format(u, v))

    if verbose:
        print('Residuals = [{}, {}]'.format(xyRes[0], xyRes[0]))
        # Print out the results
        print('origin: {}'.format(origin))
        print('u: {}\nv: {}'.format(u, v))
        print('|u|: {0[0]:.4}, |v|: {0[1]:.4} (pixels)'.format(np.linalg.norm([u, v], axis=1)))
        print('angle UV: {:.4}'.format(np.arccos(np.dot(u / np.linalg.norm(u),
                                                        v / np.linalg.norm(v))) * 180. / np.pi))

    return origin, u, v, ab.T


def refineLattice3D(or0, u0, v0, w0, pos, fraction=(1, 1, 1), max_iter=30, refine_locally=True):
    """ Refine lattice based on measurements and initial guess at lattice vectors
    
    Warning
    -------
        Tested code but not production yet. This is a work un progress.
    
    Parameters
    ----------
        or0 : 3-element tuple 
            Origin.
        u0 : tuple
            Initial u vector.    
        v0 : tuple 
            Initial v vector
        w0 : tuple 
            Initial w vector
        pos : array (number, position) (M,3)
            The set of positions to fit the lattice to.
        refine_locally : bool
            Refine locally near the origin before using all positions
            Locally is considered 2 time the larger vector of u0 or v0.
        fraction : tuple
            Site fraction to take into account peak positions inside the unit cell. FCC imaged along [100] for
            example would need site_fraction = (2, 2).
        max_iter : int
            Maximum number of iterations.
    Returns
    -------
        : tuple of 4 np.ndarray (origin, u, v, ab)
            Tuple of refined origin, u, v optimized values as np.ndarray and the
            site positions in fractions of u and v.
    
    """
    xRes0 = 1e6
    yRes0 = 1e6
    zRes0 = 1e6
    for ii in range(max_iter):

        # compute a,b values
        uv = np.array([u0, v0, w0]).T
        (ab, abRes, abX, abY) = np.linalg.lstsq(uv, (pos - or0).T, rcond=None)

        # Refine lattice
        A = np.ones([ab.shape[1], 4])
        A[:, 1] = np.round(ab[0, :])
        A[:, 2] = np.round(ab[1, :])
        A[:, 3] = np.round(ab[2, :])
        (xy_beta, xyRes, rXY, sXY) = np.linalg.lstsq(A, pos, rcond=None)

        if xyRes[0] < xRes0 and xyRes[1] < yRes0 and xyRes[2] < zRes0:
            _ = xyRes[0]
            _ = xyRes[1]
            _ = xyRes[2]
        else:
            print('----Converged----')
            break

        print('Residuals = [{0[0]}, {0[1]}, {0[2]}]'.format(xyRes))

        origin = xy_beta[0, :]
        u = xy_beta[1, :]
        v = xy_beta[2, :]
        w = xy_beta[3, :]

        return origin, u, v, w


def generateLatticeFromRefinement(origin, u, v, ab, fraction=(1, 1)):
    """ Generate lattice positions from refined output
    of refineLattice2D. This can be used to generate
    all positions when fraction is used in refineLattice2D
    
    Parameters
    ----------
        origin : tuple
            Origin. 2-tuple
        u : tuple
            Initial u vector. 2-tuple
        v : tuple
            Initial v vector. 2-tuple
        fraction : tuple
            The fraction used in refineLattice2D. 2-tuple
        ab : np.ndarray
            The set of positions in terms of u and v lattice vectors. (number, position) (M, 2)
        fraction : tuple
            The fractional coordinates for the unit cell. 2-tuple
    
    Returns
    -------
        : np.ndarray
            The lattice site positions of the given lattice vectors of shape (M, 2)
    """

    lat = np.zeros_like(ab)

    lat[:, 0] = origin[0] + np.round(ab[:, 0] * fraction[0]) * (u[0] / fraction[0]) + \
                np.round(ab[:, 1] * fraction[1]) * v[0] / fraction[1]
    lat[:, 1] = origin[1] + np.round(ab[:, 0] * fraction[0]) * (u[1] / fraction[0]) + \
                np.round(ab[:, 1] * fraction[1]) * v[1] / fraction[1]

    return lat


def lattice2D_2(u, v, a, b, xy0, numPoints):
    """ A modified version of lattice2D which uses non-normalized u, v vectors

    """

    Y2D, X2D = np.meshgrid(range(numPoints[1]), range(numPoints[0]), indexing=default_indexing)

    X = X2D.ravel()
    Y = Y2D.ravel()

    xy = np.zeros((np.prod(numPoints), 2))
    xy[:, 0] = xy0[0] + a * X * u[0] + b * Y * v[0]
    xy[:, 1] = xy0[1] + a * X * u[1] + b * Y * v[1]

    return xy


def lattice3D_2(u, v, w, a, b, c, xyz0, numPoints):
    """ A modified version of lattice3D which uses non-normalized u,v vectors

    """

    Z3D, Y3D, X3D = np.meshgrid(range(numPoints[1]), range(numPoints[0]), range(numPoints[0]),
                                indexing=default_indexing)

    X = X3D.ravel()
    Y = Y3D.ravel()
    Z = Z3D.ravel()

    xyz = np.zeros((np.prod(numPoints), 3))
    xyz[:, 0] = xyz0[0] + a * X * u[0] + b * Y * v[0] + c * Z * w[0]
    xyz[:, 1] = xyz0[1] + a * X * u[1] + b * Y * v[1] + c * Z * w[1]
    xyz[:, 2] = xyz0[1] + a * X * u[2] + b * Y * v[2] + c * Z * w[2]

    return xyz


def fit_peaks_gauss2d(image, peaks, cutOut, init, bounds, remove_edge_peaks=True):
    """ Fit peaks to a 2D Gaussian. The Gaussian function is gaussND.gauss2D()
    
    Todo: Write tests. This function is copied from test_peakFind.ipynb
    Todo: Add example test_peakFind.ipynb from my jupyter notebook folder
    
    Parameters
    ----------
        image : np.ndarray
            The 2D image with the intensities to fit to
            
        peaks : np.ndarray
            The peaks in a ndarray of shape (N,2) where N is the number
            of peaks. Should match output of peakFind.peakFind2D.
            
        cutOut : int
            The size (+/-) of the region around each peak to fit.
            
        init : tuple
            A tuple of initial sigma values in x and y. (sigma_x, sigma_y)
        
        bounds : tuple
            A tuple of 2 tuples indicating the upper and lower bounds for the fitting.
            See scipy.optimize.curve_fit for more information. Ordered as
            ( (center_x_low, center_y_low, sigma_x_low, sigma_y_low), 
            (center_x_high, center_y_high, sigma_x_high, sigma_y_high))
        
        remove_edge_peaks : bool, default = True
            Peak positions at the edge of the image are set to NaN and are removed.
            If you want to have those positions returned as nan set this to False.
        
    Note
    ----
        This function uses np.meshgrid and indexing='ij' internally.
    
    Returns
    -------
        : tuple
            Returns a tuple of 3 arrays:
                [0] optimized peak positions as a (M, 2) ndarray. M <= N. Peaks on the edge of the image are removed.
                [1] peak intensity value interpolated at the optimized peak position.
                [2] the fitting values for each peak
    """

    # Indexes that sort intensity high to low
    peaks_sort_values = np.argsort(image[peaks[:, 0], peaks[:, 1]])

    # X,Y,sig_x,sig_y
    # NAN means peak near an edge and can not be fit
    optPointsNaN = np.zeros((peaks_sort_values.shape[0], 4))
    fittingValues = np.zeros_like(optPointsNaN)

    optINaN = np.zeros((optPointsNaN.shape[0]))  # holds the optimized intensity
    optI_meanNaN = np.zeros_like(optINaN)

    # Setup cutout coordinates around each peak
    X2D, Y2D = np.meshgrid(np.arange(-cutOut, cutOut + 1, 1),
                           np.arange(-cutOut, cutOut + 1, 1), indexing=default_indexing)

    for ii, index in enumerate(peaks_sort_values):

        curX = int(peaks[index, 0])
        curY = int(peaks[index, 1])

        # Check if current point + fit region is within the bounds of the volume
        if (curX >= cutOut) & (curY >= cutOut) & \
                (curX < image.shape[0] - cutOut) & \
                (curY < image.shape[1] - cutOut):

            curIm = np.float32(image[curX - cutOut:curX + cutOut + 1,
                               curY - cutOut:curY + cutOut + 1])
            curIm_norm = curIm - curIm.min()  # subtract the min
            curIm_norm = curIm_norm / curIm_norm.max()  # normalize to 1

            # Calculate error for each point based on intensity
            y_error = 1 / np.sqrt(curIm_norm.ravel() + 0.00001)  # add a small offset to remove divide by zero
            
            optP2D, optCov2D = opt.curve_fit(gaussND.gauss2D_FIT, (X2D, Y2D), curIm_norm.ravel(),
                                             p0=(0, 0, init[0], init[1]), bounds=bounds,
                                             sigma=y_error)

            fittingValues[ii, :] = optP2D

            optPointsNaN[ii, :] = np.array([curX + optP2D[0], curY + optP2D[1],
                                            optP2D[2], optP2D[3]])  # X,Y,sigma_x,sigma_y
        else:
            # Peak is near an edge
            optPointsNaN[ii, :] = np.array([np.nan] * 4)
            fittingValues[ii, :] = np.array([np.nan] * 4)
            optINaN[ii] = np.nan
            optI_meanNaN[ii] = np.nan

    # Get the intensity of each optimize peak position
    optINaN = ndimage.map_coordinates(image.astype('f'),
                                      ((optPointsNaN[:, 0]),
                                       (optPointsNaN[:, 1])),
                                      cval=np.nan)

    # Remove NANs where peak is at the edge if desired
    if remove_edge_peaks:
        optPoints = optPointsNaN[~np.isnan(optPointsNaN[:, 0]), :]
        optI = optINaN[~np.isnan(optINaN)]
        fittingValues = fittingValues[~np.isnan(optINaN)]
    else:
        optPoints = optPointsNaN
        optI = optINaN
        fittingValues = fittingValues

    return optPoints, optI, fittingValues


def fit_peaks_gauss3d(volume, peaks, cutOut, init, bounds, remove_edge_peaks=True):
    """ Fit peaks to a 2D Gaussian. The Gaussian function is gaussND.gauss2D()

    Todo: Write tests. This function is copied from test_peakFind_3d.ipynb

    Parameters
    ----------
        volume : np.ndarray
            The 3D image with the intensities to fit to.

        peaks : np.ndarray
            The peaks in a ndarray of shape (N, 3) where N is the number
            of peaks. Should match output of peakFind.peakFind3D.

        cutOut : int
            The size (+/-) of the region around each peak to fit.

        init : tuple
            A 3-tuple of initial sigma values in x and y. (sigma_x, sigma_y, sigma_z)

        bounds : tuple
            A 2-tuple of 6-tuples indicating the upper and lower bounds for the fitting.
            See scipy.optimize.curve_fit for more information. Ordered as
            ( (center_x_low, center_y_low, center_z_low, sigma_x_low, sigma_y_low, sigma_z_low),
            (center_x_high, center_y_high, center_z_high, sigma_x_high, sigma_y_high, sigma_z_high))

        remove_edge_peaks : bool, default = True
            Peak positions at the edge of the image are set to NaN and are removed.
            If you want to have those positions returned as nan set this to False.

    Note
    ----
        This function uses np.meshgrid and indexing='ij' internally.

    Returns
    -------
        : tuple
            Returns a tuple of 3 arrays:
                [0] optimized peak positions as a (M, 3) ndarray. M <= N. Peaks on the edge of the image are removed
                unless otherwise specified. Then the edge peaks are returned as NAN values.
                [1] peak intensity value interpolated at the optimized peak position.
                [2] the local fitting values for each peak.
    """
    # Sort from highest to lowest intensity
    validPeaks_sort_values = np.argsort(volume[peaks[:, 0], peaks[:, 1], peaks[:, 2]])

    # X,Y,Z,sig_x,sig_y,sig_z
    # NAN means peak near an edge and can not be fit
    optPointsNaN = np.zeros((len(validPeaks_sort_values), 6))
    fittingValues = np.zeros_like(optPointsNaN)

    # Set up fitting static variables
    local_region = np.arange(-cutOut, cutOut + 1, 1)
    Y3D, X3D, Z3D = np.meshgrid(local_region, local_region, local_region, indexing=default_indexing)
    bad_point = (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    for ii, index in enumerate(validPeaks_sort_values):
        curX = int(peaks[index, 0])
        curY = int(peaks[index, 1])
        curZ = int(peaks[index, 2])

        # Check if current point + fit region is within the bounds of the volume
        if (curX > cutOut + 1) & (curY > cutOut + 1) & (curZ > cutOut + 1) & (curX < volume.shape[0] - cutOut) \
                & (curY < volume.shape[1] - cutOut) & (curZ < volume.shape[2] - cutOut):

            curVol = np.float32(volume[curX - cutOut:curX + cutOut + 1,
                                curY - cutOut:curY + cutOut + 1,
                                curZ - cutOut:curZ + cutOut + 1])
            curVol_norm = curVol - curVol.min()  # subtract the min
            curVol_norm = curVol_norm / curVol_norm.max()  # normalize maximum to 1

            # Fit to a 3D Gaussian to the area around the peak
            y_error = 1 / np.sqrt(curVol_norm.ravel() + 0.00001)  # add error
            optP3D, optCov3D = opt.curve_fit(gaussND.gauss3D_FIT, (X3D, Y3D, Z3D), curVol_norm.ravel(),
                                             p0=(0, 0, 0, *init), bounds=bounds, sigma=y_error)
            fittingValues[ii, :] = optP3D

            # Enter values into output array as X,Y,Z,sigma_x,sigma_y,sigma_z
            optPointsNaN[ii, :] = np.array((float(curX) + optP3D[0], float(curY) + optP3D[1], float(curZ) + optP3D[2],
                                            optP3D[3], optP3D[4], optP3D[5]))

        else:
            # if the peak is near an edge then ignore it and set its values to NAN
            optPointsNaN[ii, :] = np.array(bad_point)
            fittingValues[ii, :] = np.array(bad_point)

    # Get the intensity of each optimize peak position
    optINaN = ndimage.map_coordinates(volume.astype('f'),
                                      ((optPointsNaN[:, 0]),
                                       (optPointsNaN[:, 1]),
                                       (optPointsNaN[:, 2])),
                                      cval=np.nan)

    # Remove NANs where peak is at the edge if desired
    if remove_edge_peaks:
        optPoints = optPointsNaN[~np.isnan(optPointsNaN[:, 0]), :]
        optI = optINaN[~np.isnan(optINaN)]
        fittingValues = fittingValues[~np.isnan(optINaN)]
    else:
        optPoints = optPointsNaN
        optI = optINaN
        fittingValues = fittingValues

    return optPoints, optI, fittingValues


def match_lattice_peaks(peaks, u, v, origin):
    """ Find the matching lattice points to the experimental peaks. This is useful to generate all fitted lattice points
    for a set of experimental peaks. This can then be used directly to calculate displacements for example.

    Parameters
    ----------
    peaks : ndarray
        The experimental peaks as a (num_peaks, 2) ndarray
    u, v : tuple
        The u and v vectors for the fitted lattice
    origin : tuple
        The origin of the fitted lattice for u, v vectors

    Returns
    -------
    : ndarray
        An array of shape (num_peaks, 2) where each coordinate is the closest lattice point for each input peak.
    """

    uv = np.asarray((u, v))
    ab_nearest = np.dot(peaks - origin, np.linalg.inv(uv))

    return ab_nearest


def latticeDisplacements(peaks, u, v, origin):
    """ Find the displacements of the experimental peaks to the fitted lattice parameters

    Parameters
    ----------
    peaks : ndarray
        The experimental peaks with shape (num_peaks, 2). For example, the output of fit_peaks_gauss2d.
    u, v : tuple
        The u and v vectors for the fitted lattice
    origin : tuple
        The origin of the fitted lattice for u, v vectors

    Returns
    -------
    : ndarray
        The displacement of each peak from the expected lattice position.
    """
    ab_nearest = match_lattice_peaks(peaks, u, v, origin)

    uv = np.asarray((u, v))

    p0_disp_linalg = np.zeros(peaks.shape)
    # p0_disp_r_theta_linalg = np.zeros(peaks.shape)
    for ii in range(peaks.shape[0]):
        pp = peaks[ii, :]
        rref = np.dot(np.round(ab_nearest[ii, :]), uv) + origin
        p0_disp_linalg[ii, :] = pp - rref
        # p0_disp_r_theta_linalg[ii, :] = (np.sqrt((pp[0] - rref[0]) ** 2 + (pp[1] - rref[1]) ** 2),
        #                                  np.arctan2((pp[1] - rref[1]), (pp[0] - rref[0])))
    return p0_disp_linalg


def calculate_unit_cell(image, lattice, u, v, unit_cell_size):
    """ IN PROGRESS
    Calculate a unit cell using the input lattice and lattice parameters. Any unit cell at the edge of the image is not
    used as part of the unit cell will be invalid.

    Parameters
    ----------
    image : ndarray
        The image containing the unit cell in a regular pattern (i.e. a lattice)
    lattice: ndarray
        The starting position of each unit cell in the image.
    u, v : tuple
        The u and v vectors for the fitted lattice
    unit_cell_size : int or tuple
        If int then this is the number of pixels in the unit cell along the u vector. The size of the v vector will be
        automatically calculated to make the unit cell pixels square. If tuple then the number of pixels along each
        u and v dimension are used directly.

    Returns
    -------
    : ndarray
        An array of size unit_cell_size which contains the average value of each position in the unit cell for each
        peak position and u,v, lattice point.
    """
    print('Warning, untested code. Work in progress.')
    if isinstance(unit_cell_size, int):
        num_u = unit_cell_size
        num_v = int(num_u * (np.linalg.norm(v) / np.linalg.norm(u)))
    elif isinstance(unit_cell_size, tuple):
        num_u, num_v = unit_cell_size
    else:
        raise TypeError('unit_cell_size must be int or 2-tuple')

    # Ensure image is floating point for interpolation
    if np.issubdtype(image.dtype, np.floating):
        image2 = image.copy()
    else:
        image2 = image.astype(np.float32)

    uu = [ii / num_u for ii in u]
    vv = [ii / num_v for ii in v]

    # Create a centered set of sub-lattice coordinates
    sub_lattice = lattice2D(uu, vv, 1, 1, (0, 0), (num_u, num_v))  # starts at (0, 0). Then offset for each peak.

    unit_cell = np.zeros((num_v, num_u), dtype=image2.dtype)
    cur_cell = np.zeros((num_v * num_u,), dtype=image2.dtype)
    num = 0  # number of unit cells used (edges are not used)
    for ii, peak in enumerate(lattice):
        cur_XX = sub_lattice[:, 0] + peak[0]
        cur_YY = sub_lattice[:, 1] + peak[1]
        ndimage.map_coordinates(image2, (cur_XX.ravel(), cur_YY.ravel()), order=3, output=cur_cell)
        if not np.any(np.isnan(cur_cell)):
            # Avoid unit cells on the edge of the image. Otherwise the mean is corrupted
            unit_cell += cur_cell.reshape((num_v, num_u))
            num += 1

    # Normalize by the number of unit cells used
    unit_cell /= num

    # Rotate to match the image
    return unit_cell.T
