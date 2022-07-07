from scipy import ndimage
import numpy as np


def _gen_line(p0, p1, num_points):
    """ Generate a set of points on a line between two points.
    Helper function for line_profile.


    Parameters
    ----------
        p0 : 2-tuple
            The x0 and y0 starting points of the line

        p1 : 2-tuple
            The x1 and y1 ending points of the line

        num_points : integer
            The number of points along the line

    Returns
    -------
        : tuple, x and y
            The x coordintes at which the line proifle was taken at and y the intensity values.

    """

    x = np.linspace(p0[0], p1[0], num_points)
    m = (p1[1] - p0[1]) / (p1[0] - p0[0])
    b = p0[1] - p0[0] * m
    y = m * x + b
    return x, y


def line_profile(im0, p0, p1, num_points, width=0, step=0.5):
    """ Use interpolation to measure a line scan across a 2D image. map_coordinates
    uses a different convention from numpy indexing. So the x and y positions need
    to be flipped.

    Parameters
    ----------
        im0 : ndarray
            The 2D array of the image
        p0 : tuple
            The start of the line scan (x0,y0) and (row,col)
        p1 : tuple
            The end of the line scan (x1,y1)
        num_points : int
            The number of points in the line scan
        width : int
            The width of the line scan to average over. Must be an integer.
        step : float
            The step size laterally for the interpolation. Must be < 1

    Returns
    -------
        : tuple
            A tuple containing the line profile and the positions where the interpolation was performed

    Note
    ----
        As an example: if width = 1 and step = 0.5 then there will be 3 line scans averaged.

    Example
    --------
        >> line, (xx, yy) = line_profile(image,(0,100),(175,100),50,step=0.5,width=5)
        >> fg, ax = plt.subplots(1,2)
        >> ax[0].plot(line,'*-')
        >> ax[1].imshow(image)
        >> ax[1].scatter(xx, yy)
    """

    im = im0.astype(np.float32)

    # The line scan from p0 to p1
    (x0, y0) = _gen_line(p0, p1, num_points)
    line0 = ndimage.map_coordinates(im, np.vstack((x0, y0)))  # reverse x and y to switch between row and column major

    # Average perpendicular line scans
    if width >= 1:
        x1 = []  # save all x positions
        y1 = []  # save all y positions
        x1.append(x0)
        y1.append(y0)
        x0 = x1
        y0 = y1

        # Find perpendicular direction to line scan for averaging
        v0 = np.matrix([p1[0] - p0[0], p1[1] - p0[1], 0])
        v0 = v0 / np.linalg.norm(v0)
        v2 = np.cross(v0, np.matrix([0, 0, 1]))
        v2 = v2 / np.linalg.norm(v2)
        v2 = v2[0]  # the perpendicular direction

        # Positive direction
        numLines = int(width / step / 2.0)
        for ii in range(0, numLines):
            p2 = np.array(p0) + v2[0:2] * step * ii
            p3 = np.array(p1) + v2[0:2] * step * ii
            (x, y) = _gen_line(p2, p3, num_points)
            line0 += ndimage.map_coordinates(im, np.vstack(
                (x, y)))  # reverse x and y to switch between row and column major
            x0.append(x)
            y0.append(y)
        # Negative direction
        for ii in range(0, numLines):
            p4 = np.array(p0) - v2[0:2] * step * ii
            p5 = np.array(p1) - v2[0:2] * step * ii
            (x, y) = _gen_line(p4, p5, num_points)
            line0 += ndimage.map_coordinates(im, np.vstack(
                (x, y)))  # reverse x and y to switch between row and column major
            x0.append(x)
            y0.append(y)
    line0 = line0 / (width / step + 1)

    return line0, (x0, y0)


if __name__ == '__main__':

    XX, YY = np.mgrid[0:100,0:100]
    RR = np.sqrt((XX-50)**2 + (YY-50)**2)

    profile, (xx, yy) = line_profile(RR, (0,0), (99, 99), 100)
    print(profile)
