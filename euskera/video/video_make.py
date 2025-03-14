
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from IPython.display import HTML

# Get the parent directory dynamically
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ""))
sys.path.append(parent_dir)

import plots as pl


########### Make a video
#############################################################################
class Visualization:

    def __init__(self, address, figData, solAnalit=None, info=False):
        if info:
            print(f"Usando dirección {address}")
        
        # Atributos de instancia
        self.address = address
        self.figData = figData
        #self.solAnalit = solAnalit

    # Method
    
    ### Globals
    def loadData(self, address=None):
        """
        Load the data from the self.address (default) or from a specify address (default None)
        """
        try:
            # Use the provided address or fall back to self.direccion
            file_path = address if address else getattr(self, "address", None)
        
            if not file_path:
                raise ValueError("No valid address provided to load the data.")
        
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"The file '{file_path}' does not exist.")
        
            array_data = np.load(file_path)
            return array_data
    
        except Exception as e:
            print(f"Error loading data: {e}")
            return None
    
    def extDataSave(self, nameP, name2=False,
                    coord=['x', 'y', 'z'], namegrid='end_grid.npz',
                    format='.npz', address=None):
        """ 
        Extracts data from the specified grid file and loads additional data.

        Parameters:
        - coord: List of coordinate labels (default: ['x', 'y', 'z']).
        - namegrid: Name of the grid file to load.
        - address: Directory where the grid file is stored.

        Returns:
        - A list containing the grid data (as a NumPy array) and additional loaded data.
        """

        file_path = address if address else getattr(self, "address", None)
        if not file_path:
            raise ValueError("No valid address provided to load the data.")
        
        address_coord = os.path.join(file_path, namegrid)

        # Load grid data
        grid_temp = self.loadData(address_coord)
        try:
            grid = np.array([grid_temp[axe] for axe in coord])
        except KeyError as e:
            raise KeyError(f"Missing expected coordinate key in grid data: {e}")

        # Load additional array data
        address_coord = os.path.join(file_path, 'end_' + nameP + format)
        array_data = self.loadData(address_coord)
        
        if name2:
            address_coord = os.path.join(file_path, 'end_' + name2 + format)
            array_data2 = self.loadData(address_coord)
            
        return [grid, array_data, array_data2] if name2 else [grid, array_data]
    

    #### 2D video
    
    def imagshow02D(self, ax, data, datl, xlim=(-1, 1), ylim=(-1, 1), 
                    cmapint=['#050505', '#f0784d'], xlimE=(-1, 1),
                    ylimE=(-1, 1), show=False, text=None):
        
        # Unpacking input data
        perf2d, t0 = data
        perf2d = perf2d.T
        
        perf, xi, _ = datl
        
        [xminE, xmaxE], [yminE, ymaxE] = xlimE, ylimE
        vmax = np.abs(perf2d).max()
        vmin = np.abs(perf2d).min()
    
        # Convert color list to colormap
        cmap, norm = pl.colorBar_and_normaliz(vmax, vmin, cmapStr=False, cmapint=cmapint)
        
        frame = ax[0].imshow(perf2d, cmap=cmap, extent=(xminE, xmaxE, yminE, ymaxE), norm=norm, origin='lower', zorder=1)
        
        perfmax = max(perf) if np.any(perf) else 1  # Avoid division by zero
        frame1, = ax[1].plot(xi, perf/perfmax, ls=' ', lw=.5, c='k', zorder=1)
        
        color = perf/perfmax
        frame2 = pl.colored_line(xi, perf/perfmax, color, ax[1], add=False, linewidth=2, cmap=cmap, zorder=10)
        ax[1].add_collection(frame2)
        
        # Compute text position
        x_pos = xlim[0] + (xlim[1] - xlim[0]) / 2
        y_pos = ylim[0] + (ylim[1] - ylim[0]) / 2 + 0.4
        tframe = ax[0].text(x_pos, y_pos, s=r'time = $%3.1f$' % t0, c='white', fontsize='small')
        
        if text:
            ax[0].text(xminE+x_pos/2, y_pos, s=text, c='white', fontsize='small')
        
            
        ax[0].set_axis_off()
        ax[0].set_xlim(xlim)
        ax[0].set_ylim(ylim)
        ax[1].set_xticks([])
        ax[1].set_yticks([])
          
        # Show plot if requested
        if show:
            plt.show() 
        
        return ax, frame, frame1, frame2, tframe, cmap, perfmax
    
    
    def frame02D(self, data, struc, xlim=None, ylim=None, show=False, Npart=None):
        """ 
        Creates the first frame of a 2D plot.

        Parameters:
        - data: Tuple containing (performance values, x-coordinates, initial time t0).
        - struc: Tuple containing (xmin, xmax, linestyle, linewidth, color).
        - xlim: Tuple specifying x-axis limits (default: None).
        - ylim: Tuple specifying y-axis limits (default: None).
        - show: Boolean flag to display the plot (default: False).

        Returns:
        - fig, ax, frame, tframe: The figure, axes, plot line, and text object.
        """
        # Unpacking input data
        try:
            perf, xi, t0 = data
            ls, lw, color = struc
        except ValueError:
            raise ValueError("Invalid input: Ensure 'data' has (perf, xi, t0) and 'struc' has (ls, lw, color).")
    
        # Ensure self.figData exists
        if not hasattr(self, 'figData'):
            raise AttributeError("self.figData is not defined. Ensure it's initialized before calling this function.")
        ###############################
    
        fig, ax = self.figData
        # perfm = np.max(perf)
        frame, = ax.plot(xi, perf, ls=ls, c=color, lw=lw)
        
        if Npart:
            ax.hlines(Npart, xmin=min(xi), xmax=max(xi), ls='-', lw=1, color='k', alpha=0.5)
            ax.text(x=min(xi)+1, y=Npart-0.5, s=r'Particle Number', fontsize='small')

        # Compute text position
        x_pos = xlim[0] + (xlim[1] - xlim[0]) / 2 if xlim else np.mean(xi) + 1
        y_pos = max(perf) - max(perf) / 8
        tframe = ax.text(x_pos, y_pos, s=r'time=$%3.2f$' % t0, fontsize='small')

        # Set axis limits
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(min(xi), max(xi))

        if ylim:
            ax.set_ylim(ylim)
        
        # Labels
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'Profile')

        # Show plot if requested
        if show:
            plt.show()

        return fig, ax, frame, tframe
    
    def updateimagshow02D(self, ind, frame, frame1, frame2, tframe,
                        ax, ti, name, name2, array_data, array_data2,
                        xi, cmap, perfmax):
        """ 
        """
        t = ti[ind]
        tframe.set_text(r'time=$%4.3f$'%t)

        perf1d = array_data[name+'%d'%ind]
        perf2d = array_data2[name2+'%d'%ind]
        perf2d = perf2d.T

        # updating axis
        frame.set_array(perf2d)
        
        frame1.set_ydata(perf1d)
        
        # perfmax = max(perf1d) if np.any(perf1d) else 1  # Avoid division by zero
        color = perf1d/perfmax
        
        ls = pl.colored_line(xi, perf1d/perfmax, color, ax[1], add=False, linewidth=2, cmap=cmap, zorder=10)
        frame2.update({'segments': ls.get_segments(),
                       #'cmap': ls.get_cmap(),
                       'array': color #ls.get_array()
                       })
        
        # updating y-lim
        #axT[1].set_ylim(min(ui2)-min(ui2)/8, max(ui2)+max(ui2)/8)
        
        return frame, frame1, frame2, tframe

    def updateframe2D(self, ind, frame, tframe,
                        ax, ti, name, array_data):
        """ 
        """
        t = ti[ind]
        #print(t)
        tframe.set_text(r'time=$%4.3f$'%t)

        ui = array_data[name+'%d'%ind]

        # updating axis
        # uimax = np.max(ui)
        frame.set_ydata(ui)

        # updating y-lim
        # ax.set_ylim(min(ui)-min(ui)/8, max(ui)+max(ui)/8)

        return frame, tframe

    def fplot2D(self, n0, grid, array_data, name, struc, xlim=(-1, 1),
                ylim=(-1, 1), xlimE=(-1, 1), ylimE=(-1, 1), interval=200,
                cmapint=['#050505', '#f0784d'], array_data2=None, name2=False,
                show=False, text=None, Npart=None):
        """
        Creates a 2D animated plot using time-dependent data.

        Parameters:
        - n0: Initial frame index.
        - grid: Grid data (assumed to be a list or NumPy array).
        - array_data: Dictionary containing time series and profiles.
        - name: Base name of profile data in `array_data`.
        - struc: Tuple containing (xmin, xmax, linestyle, linewidth, color).
        - show: Boolean flag to display the plot (default: False).

        Returns:
        - anim: Matplotlib animation object.
        """
        # Check if grid is valid
        if not isinstance(grid, (list, np.ndarray)) or len(grid) == 0:
            raise ValueError("Invalid 'grid': Expected a non-empty list or NumPy array.")
        # poner check the array_data
    
        # time-data
        ti = array_data['t']
        profile_key = f"{name}{n0}"
        prof0 = array_data[profile_key]
        data0 = [prof0, grid[0], ti[n0]]
        
        if name2:
            ti2 = array_data2['t']
            profile_key = f"{name2}{n0}"
            prof20 = array_data2[profile_key]
            data1 = [prof20, ti2[n0]]
        
        # Number of frames
        nframes = len(ti) - n0
        
        # First frame
        if name2:
            fig, ax = self.figData
            ax, frame, frame1, frame2, tframe, cmap, perfmax = self.imagshow02D(ax=ax, data=data1, datl=data0,
                                                                          xlim=xlim, ylim=ylim,
                                                                          cmapint=cmapint, xlimE=xlimE, ylimE=ylimE,
                                                                          show=show, text=text)
            # Animation
            anim = animation.FuncAnimation(
                fig, 
                self.updateimagshow02D,
                frames=range(n0, nframes),
                fargs=(frame, frame1, frame2, tframe, 
                       ax, ti2, name, name2, array_data, array_data2, grid[0], cmap, perfmax),
                interval=interval, 
                blit=False
                )
        else:
            fig, ax, frame, tframe = self.frame02D(data0, struc, xlim=xlim, ylim=ylim, show=show, Npart=Npart)
            
            # Animation
            anim = animation.FuncAnimation(
                fig, 
                self.updateframe2D,
                frames=range(n0, nframes),
                fargs=(frame, tframe, ax, ti, name, array_data),
                interval=interval, 
                blit=False
                )
        
        return anim

    ############

    #############
    def video(self, nameP, struc, nameV, name2=False, coord=['x', 'y', 'z'], address=None,
              namegrid='end_grid.npz',
              format='.npz',
              plot2D=True, plot3D=False,
              n0=0, show=False,
              xlim=(-1, 1), ylim=(-1, 1), xlimE=(-1, 1), ylimE=(-1, 1),
              cmapint=['#050505', '#f0784d'], interval=200,
              vconf=[10, 1000000, ['-vcodec', 'libx264']], save=True, text=None,  Npart=None):
        """ 
        Making a video from the data
        """
        
        # loading data
        if name2:
            grid, array_data, array_data2 = self.extDataSave(coord=coord, namegrid=namegrid,
                                                nameP=nameP, name2=name2, address=address, format=format)
        else:
            grid, array_data = self.extDataSave(coord=coord, namegrid=namegrid,
                                            nameP=nameP, name2=name2, address=address, format=format)
            array_data2 = None

        if plot2D:
            anim = self.fplot2D(n0, grid, array_data, nameP, struc,  cmapint=cmapint, array_data2=array_data2,
                                name2=name2, xlim=xlim, ylim=ylim, xlimE=xlimE, ylimE=ylimE, show=show, 
                                interval=interval, text=text, Npart=Npart)
        if plot3D:
            pass

        fps, bitrate, extra_args = vconf
        
        if save:
            anim.save(nameV+'.mp4', bitrate=bitrate, fps=fps, extra_args=extra_args)  # direc+nameV+
        else:
            animated_plot = HTML(anim.to_jshtml())
            return animated_plot
        
        return None