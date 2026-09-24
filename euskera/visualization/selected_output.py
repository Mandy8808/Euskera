"""Coordinate-aware plots for selected line and plane streams."""

import numpy as np
import matplotlib.pyplot as plt
from euskera.io.selected_output import read_output

__all__ = ["plot_output"]


def plot_output(address, name, sample=-1, component=0, ax=None):
    """Plot a stored line or plane; sample is a position, not a global index.

    For psi, display the modulus of the requested component. The returned
    Matplotlib axes can be further customized by the caller.
    """
    output = read_output(address, name)
    if output["geometry"] not in ("line", "plane"):
        raise ValueError("plot_output requires a line or plane stream")
    if not len(output["time"]):
        raise ValueError("The selected output contains no samples")
    values = output["data"][sample]
    if output["variable"] == "psi":
        values = np.abs(values[component])
    if ax is None:
        _, ax = plt.subplots()
    axes = output["axes"]
    coords = output["coordinates"]
    if len(axes) == 1:
        ax.plot(coords[axes[0]], values)
        ax.set_ylabel(f"|psi[{component}]|" if output["variable"] == "psi" else output["variable"])
    else:
        artist = ax.pcolormesh(coords[axes[0]], coords[axes[1]], values.T, shading="auto")
        ax.figure.colorbar(artist, ax=ax)
        ax.set_ylabel(axes[1])
    ax.set_xlabel(axes[0])
    ax.set_title(f"{name}: t={output['time'][sample]:.6g}, sample={output['snapshot_index'][sample]}")
    return ax
