# euskera v1.0
# plot_tools file

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
###################################################################################################


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


########### 2D-Plane
#############################################################################
def imagshow2D(data, datl=None, xlim=(-0.45, 0.45), ylim=(-0.45, 0.45),
               contour=True, imshow=True, cmapint=['#050505', '#f0784d'],
               xlimE=(-1, 1), ylimE=(-1, 1), levels=4, alpha=0.5,
               save=False, filename='rho2D.pdf', dpi=1000, out=False, v_vals=None, info=False):
    """
    Function to visualize 2D data using contour and imshow.
    """

    [xminE, xmaxE], [yminE, ymaxE] = xlimE, ylimE
    vmin, vmax = v_vals if v_vals else [np.abs(data).min(), np.abs(data).max()]
    if info:
        print('vmin, vmax', vmin, vmax)

    cmap, norm = colorBar_and_normaliz(vmax, vmin, cmapStr=False, cmapint=cmapint)

    if datl is None:
        fig, ax = plt.subplots(figsize=(8, 8))
        axes = [ax]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
        ax = axes[0]

    if contour:
        ax.contour(data, cmap=cmap, norm=norm, origin='lower', alpha=alpha,
                   extent=(xminE, xmaxE, yminE, ymaxE), levels=levels, zorder=2)

    if imshow:
        ax.imshow(data, cmap=cmap, extent=(xminE, xmaxE, yminE, ymaxE),
                  norm=norm, origin='lower', zorder=1)

    ax.set_axis_off()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if datl is not None:
        x, y = datl
        y_norm = y / np.max(y)
        axes[1].plot(x, y_norm, ls=' ', lw=.5, c='k', zorder=1)
        colored_line(x, y_norm, y_norm, axes[1], linewidth=2, cmap=cmap, zorder=10)
        axes[1].set_xticks([])
        axes[1].set_yticks([])

    if save:
        fig.savefig(filename, format='pdf', pad_inches=0.1, dpi=dpi, bbox_inches='tight')
    elif not out:
        plt.show()

    return fig, axes

def old_imagshow2D(data, datl=None, xlim=(-0.45, 0.45), ylim=(-0.45, 0.45), contour=True, imshow=True,
               cmapint=['#050505', '#f0784d'], xlimE=(-1, 1), ylimE=(-1, 1), 
               levels=4, alpha=0.5, save=False, out=False):
    """
    Function to visualize 2D data using contour and imshow.
    """
    
    [xminE, xmaxE], [yminE, ymaxE] = xlimE, ylimE
    vmax = np.abs(data).max()
    vmin = np.abs(data).min()
    
    # Convert color list to colormap
    cmap, norm = colorBar_and_normaliz(vmax, vmin, cmapStr=False, cmapint=cmapint)

    fig, axT = plt.subplots(figsize=(8, 8)) if not datl else plt.subplots(1, 2, figsize=(8, 3.5))  # (width=6, height=8)
    ax = axT if  not datl else axT[0]
    
    if contour:
        ax.contour(data, cmap=cmap, norm=norm, origin='lower', alpha=alpha, 
                   extent=(xminE, xmaxE, yminE, ymaxE), levels=levels, zorder=2)
    
    if imshow:
        ax.imshow(data, cmap=cmap, extent=(xminE, xmaxE, yminE, ymaxE), norm=norm, origin='lower', zorder=1)

    ax.set_axis_off()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if datl:
        x, y = datl
        ymax = max(y)
        axT[1].plot(x, y/ymax, ls=' ', lw=.5, c='k', zorder=1)
        color = y/ymax
        colored_line(x, y/ymax, color, axT[1], linewidth=2, cmap=cmap, zorder=10)
        axT[1].set_xticks([])
        axT[1].set_yticks([])
    
    if save:
        fig.savefig('rho2D.pdf', format='pdf', pad_inches=0.1, dpi=1000, bbox_inches='tight')
    else:
        if not out:
            plt.show()
    
    return fig, axT

########### 3D-plot
#############################################################################
def ShowPlaneProf(profData, X=None, Y=None, Z=None, indX=None, indY=None, indZ=None, save=False):
    """ 
    Function to validate and process plane profiles.
    """
    
    ####### Validate numerical parameters
    # Ensure the correct structure of the profile data
    if len(profData.shape) != 3:
        raise ValueError("The profile needs to have the structure: (resol, resol, resol)")

    # Ensure exactly two coordinate components are provided
    if sum(v is not None for v in [X, Y, Z]) != 2:
        raise ValueError("It is necessary to provide data for exactly two coordinate components.")

    # Ensure indices are either int or None, and only one index is provided
    if not all(isinstance(temp, (int, type(None))) for temp in [indX, indY, indZ]):
        raise ValueError("The indices must be integers or None.")
    if sum(v is not None for v in [indX, indY, indZ]) not in [0, 1]:
        raise ValueError("Only one index must be provided.")
    

    ####### Preparing data
    coord = np.array(["x", "y", "z"], dtype=object)
    xd, yd = [temp for temp in [X, Y, Z] if temp is not None]
    cxd, cyd = [coord[k] for k, temp in enumerate([X, Y, Z]) if temp is not None]

    # Extract profile based on provided index
    if indX is not None:
        prof = profData[indX, :, :]
    elif indY is not None:
        prof = profData[:, indY, :]
    elif indZ is not None:
        prof = profData[:, :, indZ]
    else:
        print("Showing the Z=0, X-Y plane")
        resol = profData.shape[0]
        prof = profData[:, :, resol // 2]  # Ensure integer division
        
    PlaneProf(prof, xd, yd, cxd, cyd, save=save)
    
    return None

def PlaneProf(prof, xd, yd, cxd, cyd, save=False, v_vals=None, info=False):
    """ 
    Plot a 2D proyection
    """
    
    # Generate mesh grid
    xdG, ydG = np.meshgrid(xd, yd, indexing='ij')

    # Get min/max values for normalization
    vmin, vmax = v_vals if v_vals else [np.abs(prof).min(), np.abs(prof).max()]
    if info:
        print('vmin, vmax', vmin, vmax)
    #maximo = np.max(prof)
    #minimo = np.min(prof)
    
    cmap, norm = colorBar_and_normaliz(vmax, vmin)  # Ensure function name matches
    
    ###### Plot 3D profile
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Normalize color values
    color_values = norm(prof.flatten())  # Flatten data for color mapping
    
    # Scatter plot
    ax.scatter(xdG.flatten(), ydG.flatten(), prof.flatten(), c=color_values, cmap=cmap, norm=norm, s=0.5)

    ax.set_xlabel(r'%s' % cxd)
    ax.set_ylabel(r'%s' % cyd)
    ax.set_zlabel(r'Profile')

    ax.set_xlim(np.min(xdG), np.max(xdG))
    ax.set_ylim(np.min(ydG), np.max(ydG))
    ax.set_zlim(np.min(prof), np.max(prof))
    
    # Line plots for 2D reference
    ax.plot(xd, prof[:, prof.shape[1] // 2], zs=0, zdir='y', color='k', alpha=0.8)
    ax.plot(yd, prof[prof.shape[0] // 2, :], zs=0, zdir='x', color='r', alpha=0.8)
    plt.show()

    # Save or show plot
    if save:
        fig.savefig('P0plot.pdf', format='pdf', pad_inches=0.1, dpi=1000, bbox_inches='tight')
    else:
        #ax.view_init(elev=90., azim=-90)
        plt.show()
    
    return None

def Plot3DCorrProf(X, Y, Z, profile, minpsi=1e-8, save=False):
    """ 
    Function to plot a 3D correlation profile.
    """
    
    # Generate mesh grid (Ensure correct shape)
    X, Y, Z = np.meshgrid(X, Y, Z, indexing='ij')  # Ensure correct indexing
    
    # Define max/min for normalization
    maximo = np.max(profile)
    minimo = np.min(profile)

    # Get colormap and normalization function
    cmap, norm = colorBar_and_normaliz(maximo, minimo)
    
    # Masking data based on `minpsi` threshold
    ind = profile > minpsi

    # Create figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Normalize color values
    color_values = norm(profile[ind])  # Normalize only valid data

    # Scatter plot with proper colormap
    sc = ax.scatter(X[ind], Y[ind], Z[ind], c=color_values, cmap=cmap, norm=norm, s=0.01)

    # Labels
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')

    # Set limits
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    ax.set_zlim(np.min(Z), np.max(Z))

    # Adjust view
    ax.view_init(elev=40., azim=45)

    # Save or show plot
    if save:
        fig.savefig('P0plot.pdf', format='pdf', pad_inches=0.1, dpi=1000, bbox_inches='tight')
    else:
        plt.show()

    return None

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