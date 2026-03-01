# euskera v1.0
# video make file

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
from IPython.display import HTML

import plots as pl

#############################################################################
class Visualization:

    # Constructor
    #########################################################################
    def __init__(self, address, figData, info=False):
        self.address = address
        self.fig, ax = figData
        self.ax = np.atleast_1d(ax)  # Always array-safe

        # Internal state
        self._frame_quiver = None
        self._gxi = None
        self._gyi = None
        self._step = None
        self._line_profiles = {'x': None, 'y': None, 'z': None}
        self._line_collections = {'x': None, 'y': None, 'z': None}

        if info: print(f"Using address {address}")

    # Data loading
    #########################################################################
    def loadData(self, address=None):
        """Load the data from the self.address (default) or from a specify address (default None)"""

        file_path = address if address else self.address
        if not file_path: raise ValueError("No valid address provided.")
        if not os.path.exists(file_path): raise FileNotFoundError(f"File '{file_path}' does not exist.")
        return np.load(file_path, allow_pickle=True)

    def extDataSave(self, nameL, nameS=None, coord=('x', 'y', 'z'), namegrid='end_grid.npz',
                    format='.npz', address=None):
        """Extracts data from the specified grid file and loads additional data."""

        base_path = address if address else self.address
        if not base_path: raise ValueError("No valid address provided.")

        # Load grid data
        grid_path = os.path.join(base_path, namegrid)
        grid_temp = self.loadData(grid_path)
        grid = np.array([grid_temp[axe] for axe in coord])

        # Load primary data
        array_dataL = None
        if nameL:
            data_path = os.path.join(base_path, f'end_{nameL}{format}')
            array_dataL = self.loadData(data_path)

        array_dataS = None
        if nameS:
            data_path2 = os.path.join(base_path, f'end_{nameS}{format}')
            array_dataS = self.loadData(data_path2)

        return grid, array_dataL, array_dataS

    # Polarization
    #########################################################################
    def _compute_polarization(self, psi2d):
        """Compute the density and polarization vectors from the 2D wavefunction."""

        rho_i = np.real(psi2d * np.conj(psi2d))
        density = np.sum(rho_i, axis=0)

        threshold = 0.001 * density.max()  #  # 0.1% of maximum density
        mask = density > threshold

        inv_sqrt = np.zeros_like(density)
        inv_sqrt[mask] = 1.0 / np.sqrt(density[mask])

        Ux = np.real(psi2d[0]) * inv_sqrt
        Uy = np.real(psi2d[1]) * inv_sqrt
        return density, Ux, Uy

    # First frame creation
    #########################################################################
    def imagshow02D(self, dataS, dataL, xlim, ylim, cmapint, plot2D):
        """Create the initial frame for the 2D performance and 1D profile."""

        panel = 0
        # -------- 2D surface --------
        frame = None
        if dataS:
            perfPsi2d, t0, xi, yi = dataS

            if plot2D == "polarization":
                density, Ux, Uy = self._compute_polarization(perfPsi2d)
                perf2d = density
                self._gxi, self._gyi = np.meshgrid(xi, yi, indexing='ij')
            else:
                perf2d = perfPsi2d
            perf2d = perf2d.T

            vmax = np.abs(perf2d).max()
            vmin = np.abs(perf2d).min()
            cmap, norm = pl.colorBar_and_normaliz(vmax, vmin, cmapStr=False, cmapint=cmapint)

            frame = self.ax[panel].imshow(perf2d, cmap=cmap, 
                                 extent=(xi.min(), xi.max(), yi.min(), yi.max()),
                                 norm=norm, origin='lower'
                                 )

            # -------- quiver --------
            if plot2D == "polarization":
                grid_size = len(xi)
                self._step = max(2, grid_size // 200)  # Adjust step based on grid size, with a minimum of 2
                self._frame_quiver = self.ax[panel].quiver(
                    self._gxi[::self._step, ::self._step],
                    self._gyi[::self._step, ::self._step],
                    Ux[::self._step, ::self._step],
                    Uy[::self._step, ::self._step],
                    color="white",
                    scale=0.8,
                    pivot='middle',
                    scale_units='xy',
                    width=0.005
                )

            self.ax[panel].set_axis_off()
            self.ax[panel].set_xlim(xlim)
            self.ax[panel].set_ylim(ylim)
            panel += 1

        # -------- 1D profile --------
        perfmax = None
        if dataL:
            perf1d, xi1d, t0 = dataL
            if perf1d.shape != xi1d.shape:
                print(f"Warning: perf1d shape {perf1d.shape} does not match xi1d shape {xi1d.shape}. You introduced the Psi.")

            perfmax = np.abs(perf1d).max()
            if perfmax == 0: perfmax = 1.0  # Avoid division by zero
            ydata = np.abs(perf1d) / perfmax

            # Set up the color mapping for the 1D profile if needed
            if not dataS:  # Only set color mapping if 2D data is not present
                vmax = np.abs(ydata).max()
                vmin = np.abs(ydata).min()
                cmap, norm = pl.colorBar_and_normaliz(vmax, vmin, cmapStr=False, cmapint=cmapint)

            ydata = np.atleast_2d(ydata)  # Ensure ydata is array-like
            for axe, data in zip(['x', 'y', 'z'], ydata):
                self._line_profiles[axe], = self.ax[panel].plot(xi1d, data, lw=0.5)

                # Create LineCollection for the 1D profile
                points = np.array([xi1d, data]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                self._line_collections[axe] = LineCollection(segments, cmap=cmap, norm=norm)
                self._line_collections[axe].set_array(data)
                self.ax[panel].add_collection(self._line_collections[axe])

            self.ax[panel].set_axis_off()
           
        if not dataS and not dataL: raise ValueError("At least one of dataS or dataL must be provided.")

        # Add time annotation
        panel = 0  # Put time annotation on the first panel
        xmin, xmax = self.ax[panel].get_xlim()
        ymin, ymax = self.ax[panel].get_ylim()

        x_pos = xmin + 0.1 * (xmax - xmin)  # Position to the left of the plot 
        y_pos = ymax - 0.1 * (ymax - ymin)  # Adjust y_pos to be above the top of the plot
        tframe = self.ax[panel].text(x_pos, y_pos, f"time = {t0:.3f}", color='red', fontsize=8, verticalalignment='top')

        return frame, tframe, perfmax
    
    # Update function for animation
    #########################################################################
    def updateimagshow02D(self, ind, frame, tframe, ti, nameL, nameS,
                          array_dataL, array_dataS, xi, perfmax, plot2D):
        """Update the 2D performance and 1D profile for the animation."""

        frame, tframe = frame, tframe

        if tframe: tframe.set_text(f"time = {ti[ind]:.3f}")

        # Update 2D surface
        panel = 0
        if frame:
            if plot2D == "polarization":
                psi2d = array_dataS[f"{nameS}{ind}"]
                perf2d, Ux, Uy = self._compute_polarization(psi2d)

                self._frame_quiver.set_UVC(
                    Ux[::self._step, ::self._step],
                    Uy[::self._step, ::self._step]
                )
            else:
                if array_dataS: perf2d = array_dataS[f"{nameS}{ind}"]
            frame.set_array(perf2d.T)
            panel += 1

        # Update 1D line
        if array_dataL:
            perf1d = array_dataL[f"{nameL}{ind}"]
            ydata = np.abs(perf1d) / perfmax

            ydata = np.atleast_2d(ydata)
            for axe, data in zip(['x', 'y', 'z'], ydata):
                self._line_profiles[axe].set_ydata(data)

                # Update LineCollection
                points = np.array([xi, data]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                self._line_collections[axe].set_segments(segments)
                self._line_collections[axe].set_array(data)

            # Dynamic ylim (only ground)
            ymin = ydata.min()
            ymax = ydata.max()

            if ymin != ymax:
                margin = 0.05 * (ymax - ymin)
                new_ymin = ymin  # - margin
                new_ymax = ymax + margin
                current_ymin, current_ymax = self.ax[panel].get_ylim()

                self.ax[panel].set_ylim(
                    min(new_ymin, current_ymin),
                    max(new_ymax, current_ymax)
                    )

        return frame, tframe

    # Animation builder
    #########################################################################
    def fplot2D(self, n0, grid, array_dataL, nameL, plot2D,
                xlim, ylim, cmapint, interval, array_dataS, nameS, dt):
        """Build the animation for the 2D performance and 1D profile."""

        dataL = None
        if nameL:
            ti = array_dataL['t'] * dt
            prof0 = array_dataL[f"{nameL}{n0}"]
            dataL = [prof0, grid[0], ti[n0]]

        dataS = None
        if nameS:
            ti2 = array_dataS['t'] * dt
            prof20 = array_dataS[f"{nameS}{n0}"]
            dataS = [prof20, ti2[n0], grid[0], grid[1]]

        time = ti2 if dataS else ti  # Use the time array from 2D data if available, otherwise use 1D time array
        frame, tframe, perfmax = self.imagshow02D(dataS, dataL, xlim, ylim, cmapint, plot2D)
        anim = animation.FuncAnimation(self.fig, self.updateimagshow02D, frames=range(n0, len(time)),
                                       fargs=(frame, tframe, time,
                                              nameL, nameS, array_dataL,
                                              array_dataS, grid[0], perfmax, plot2D),
                                              interval=interval, blit=False)
        return anim

    # Main video creation function
    #########################################################################
    def video(self, nameL, nameS, save_name, coord=('x', 'y', 'z'),
              address=None, plot2D="polarization", n0=0,
              xlim=(-1, 1), ylim=(-1, 1), cmapint=('#050505', '#f0784d'),
              interval=50, vconf=(20, 2000000, ['-vcodec', 'libx264']),
              save=True, dt=1.0):
        """Main function to create the video from the data."""

        grid, array_dataL, array_dataS = self.extDataSave( nameL=nameL, nameS=nameS, coord=coord, address=address)

        anim = self.fplot2D(n0, grid, array_dataL, nameL, plot2D, xlim, ylim, cmapint, interval, array_dataS, nameS, dt)

        if save:
            fps, bitrate, extra_args = vconf
            anim.save(f"{save_name}.mp4", fps=fps, bitrate=bitrate, extra_args=extra_args)
            return None
        
        html_video = anim.to_jshtml()
        return HTML(f"""
                    <div style="width:600px; margin:auto;">
                    {html_video}
                    </div>
                    """)