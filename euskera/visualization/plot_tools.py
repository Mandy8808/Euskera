"""Reusable plotting utilities shared by simulation visualizations."""

import warnings
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection

########### Color bar
#############################################################################
def colorBar_and_normaliz(maximo, minimo, cmapStr=False, cmapint=None):
    """
    Creates a color map and normalization scale for visualizations.
    """

    # Ensure valid input ranges
    if minimo >= maximo:
        raise ValueError("minimo must be less than maximo")

    # Creating normalized color bar
    if cmapint:
        cmap = LinearSegmentedColormap.from_list('', cmapint, N=2000)
    elif cmapStr:
        cmap = LinearSegmentedColormap.from_list('',
                                             [(0, '#c23838'),
                                              (0.25, '#f0784d'),
                                              (0.5, '#337512'),
                                              (0.75, '#1dc4ab'),
                                              (1, '#2681ab')], N=2000)
    else:
        cmap = LinearSegmentedColormap.from_list('', ['#f0784d', '#2681ab'])

    norm = mcolors.Normalize(vmin=minimo, vmax=maximo)

    return cmap, norm

########### Line with gradient color
#############################################################################
# https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
def colored_line(x, y, c, ax, add=True,**lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    # Validate inputs
    x = np.asarray(x)
    y = np.asarray(y)
    c = np.asarray(c)

    if len(x) != len(y) or len(x) != len(c):
        raise ValueError("x, y, and c must all have the same length.")

    # Check for overridden 'array' in kwargs
    if "array" in lc_kwargs:
       print('WARNING: The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    # Create and add LineCollection
    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    if add:
        return ax.add_collection(lc)
    else:
        return lc
