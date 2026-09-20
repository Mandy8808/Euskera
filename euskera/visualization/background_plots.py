import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from .configuration import get_colors
###################################
#### Illustrative plot structure
###################################
def plotUsingPerf(ax, perf, ls='-',
                  xscale='linear', yscale='linear',
                  color=None):
    """
    Plot performance profiles on given axes.

    Parameters:
    ax : matplotlib.axes.Axes or array-like of Axes
        Axes object(s) to plot on.
    perf : list or 2D array
        First element (perf[0]) is assumed to be the radial values.
        Remaining elements (perf[1:], indexed by k) are plotted against perf[0].
    """

    rad = perf[0]
    if color is None:
        color = get_colors()[0]

    # Ensure ax is iterable
    if not isinstance(ax, (list, np.ndarray)):
        ax = [ax]  # Convert single AxesSubplot into a list

    # Check if perf has enough data
    if len(perf) - 1 != len(ax):
        raise ValueError(f"Mismatch: {len(perf)-1} data series but {len(ax)} axes.")

    # Plot each performance metric
    for k, axi in enumerate(ax, start=1):
        axi.plot(rad, perf[k], ls=ls, c=color, label=f'Profile {k}')

        axi.set_xlabel(r"$r$")  # Adjust label as needed
        #axi.set_ylabel(f'Profile {k}')

        axi.legend(frameon=False)
        axi.set_xscale(xscale)
        axi.set_yscale(yscale)

    return ax

def plotUsingDiscSol(perf1, perf2=None, perf3=None, xscale='linear', yscale='linear'):
    """
    Plot the solution using discrete solutions.

    Note: The first position of perf1 (perf1[0]) is assumed to contain the radii data.

    Note: scale: 'linear', 'log', 'symlog', 'logit'
    """
    col = get_colors()

    ncol = len(perf1) - 1  # Adjusting for radii data

    if ncol <= 0:  # Safety check
        print("Error: Solution data is empty or invalid.")
        return

    _, ax = plt.subplots(ncols=ncol, figsize=(4.5 * ncol, 4.5), gridspec_kw=dict(hspace=0.0, wspace=.28))

    # Ensure ax is a list if there's only one subplot
    if ncol == 1:
        ax = [ax]

    plotUsingPerf(ax, perf1, ls='-', color=col[0], xscale=xscale, yscale=yscale)
    if perf2 is not None:
        plotUsingPerf(ax, perf2, ls='--', color=col[1], xscale=xscale, yscale=yscale)
    if perf3 is not None:
        plotUsingPerf(ax, perf3, ls=':', color=col[2], xscale=xscale, yscale=yscale)

    plt.show()
    return _, ax

def plotUsingSol(sol):
    """
    Plot the solution using a solve_ivp object.
    """
    numFields = len(sol.y)//4  # Integer division

    if numFields == 0 or sol.t.size == 0:  # Safety check
        print("Error: Solution data is empty or invalid.")
        return

    _, ax = plt.subplots(figsize=(6, 4.5))

    # Loop through each field
    y0val = []
    for k in range(0, 2 * numFields, 2):
        ax.plot(sol.t, sol.y[k], label=fr'$\sigma_{{{k//2 + 1}}}$')
        y0val.append(sol.y[k][0])

    # ymax and ymin calculations
    ymax = np.max(y0val) if y0val else 0
    ymin = np.min([np.min(sol.y[k]) for k in range(0, 2 * numFields, 2)]) if numFields > 0 else 0

    ax.set_xlim(0, sol.t[-1])
    ax.set_ylim(ymin - 0.02, ymax + 0.02)
    ax.legend(frameon=False)
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\sigma$")

    plt.show()
    return _, ax

def plotPerf(U0, rTmax, arg, px=True, py=True,
             pz=True, lim=True, Rtol=1e-9, Atol=1e-10, met='RK45'):
    """
    Plot performance metrics based on numerical integration results.

    Parameters
    ----------
    U0 : array-like
        Initial conditions for the system.
    rTmax : float
        Maximum value of the independent variable `r`.
    LambT : float
        Parameter for the system dynamics.
    px, py, pz : bool, optional
        Flags to plot results for x, y, z dimensions. Defaults are True.
    lim : bool, optional
        Whether to set dynamic axis limits. Default is True.

    Returns
    -------
    None
    """
    # Import lazily: backgrounds imports this module for plotting helpers.
    from euskera.backgrounds.systems import systemMultifrequency

    # Numerical integration setup
    rmin, rmax = 0, rTmax
    Nptos = 500
    rspan = np.linspace(rmin, rmax, Nptos)

    try:
        solution = solve_ivp(systemMultifrequency, [rmin, rmax], U0, t_eval=rspan,
                             args=(arg,), method=met, rtol=Rtol, atol=Atol)
    except Exception as e:
        print(f"Error during integration: {e}")
        return

    # Choosing profiles
    numFields = len(solution.y)//4  # Integer division
    if numFields == 1:
        warnings.warn("The py, pz ccomponents was off because the number of fields is %d"%numFields)
        py, pz = False, False
    elif numFields == 2:
        warnings.warn("The pz component was off because the number of fields is %d"%numFields)
        pz = False

    # Extract solution components
    t = solution.t
    y = solution.y

    # Plot setup
    _, ax = plt.subplots(figsize=(6, 4.5))

    # Plot dimensions based on flags
    if px:
        ax.plot(t, y[0], color='#1f77b4', label=fr'$\sigma_x^{(0)} = {y[0][0]:.2f}$')
    if py:
        ax.plot(t, y[2], color='#ff7f0e', label=fr'$\sigma_y^{(0)} = {y[2][0]:.2f}$')
    if pz:
        ax.plot(t, y[4], color='k', label=fr'$\sigma_z^{(0)} = {y[4][0]:.2f}$')

    ax.legend(frameon=False)

    # Dynamic axis limits
    if lim:
        plotted = [y[index] for index, enabled in ((0, px), (2, py), (4, pz))
                   if enabled and index < len(y)]
        if plotted:
            ymin = np.min(plotted)
            ymax = np.max(plotted)
            margin = 0.1 * max(abs(ymin), abs(ymax), 1e-12)
            ax.set_ylim(ymin - margin, ymax + margin)
        ax.set_xlim(0, t[-1])

    # Additional plot settings
    ax.hlines(y=0, xmin=0, xmax=t[-1], linestyle='--', linewidth=0.5, color='k')
    ax.set_ylabel(r'$\sigma_i^{(0)}$')
    ax.set_xlabel(r'$r$')
    # ax.set_yscale('log')

    # Display the plot
    plt.show()
    return _, ax

