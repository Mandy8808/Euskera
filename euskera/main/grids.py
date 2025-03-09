# euskera v1.0
# grids file

import numpy as np
import numexpr as ne
#############################################################################


########### Functions that generan the real/complex grids
#############################################################################
def meshgrid(*xi):
    ndim = len(xi)
    s0 = (1,) * ndim
    output = [np.asanyarray(x).reshape(s0[:i] + (-1,) + s0[i + 1:]) for i, x in enumerate(xi)]
    return output

def RealGrid(gridlength=10, resol=128):
    """
    Generate a 3D grid and compute radial distances from the origin.

    Parameters:
        gridlength (float): The physical length of the grid in one dimension. Default 10
        resol (int): The number of points in each dimension (resolution). Default 128

    Returns:
        tuple: A tuple containing:
            - A list of 3 arrays (xarray, yarray, zarray) representing the grid coordinates.
            - An array (distarray) containing the radial distances from the origin for each point.
    """
    # Generate a vector of grid points centered around zero
    gridvec = np.linspace(
        -gridlength/2.0 + gridlength/(2.0*resol),
         gridlength/2.0 - gridlength/(2.0*resol),
         resol
    ) 
    
    # Generate a sparse 3D grid (saves memory for large grids)
    xarray, yarray, zarray = meshgrid(gridvec, gridvec, gridvec)  # a little more fast
                            # np.meshgrid(gridvec, gridvec, gridvec, sparse=True, indexing='ij')
    
    # Compute radial distances using numexpr for performance
    distarray = ne.evaluate("(xarray**2 + yarray**2 + zarray**2)**0.5")
    
    return [xarray, yarray, zarray], distarray


def KGrid(gridlength=10, resol=128, realspace=False):
    """
    Generate k-space coordinates and compute their squared magnitudes.

    Parameters:
        gridlength (float): Physical length of the grid in one dimension.
        resol (int): Resolution of the grid (number of points in each dimension).
        realspace (boold -> default False): Reserved for future functionality (currently unused).

    Returns:
        tuple: A tuple containing:
            - A list of 3 arrays (kxarray, kyarray, kzarray) representing the k-space coordinates.
            - An array (karray2) containing the squared magnitudes of the k-space vectors.
    """
    # Calculate the spacing in real space
    razonNum = float(gridlength/resol)
    
    # Generate the k-space vector using FFT frequencies
    kvec = 2 * np.pi * np.fft.fftfreq(resol, d=razonNum)
    krealvec = 2 * np.pi * np.fft.rfftfreq(resol, d=razonNum) if realspace else None 
    
    # Create sparse 3D arrays for k-space coordinates
    kxarray, kyarray, kzarray = meshgrid(kvec, kvec, kvec) if not realspace else meshgrid(kvec, kvec, krealvec)
                                #np.meshgrid(kvec, kvec, kvec, sparse=True, indexing='ij') if not realspace else\
                                #            np.meshgrid(kvec, kvec, krealvec, sparse=True, indexing='ij')
    
    # Compute squared magnitudes of the k-space vectors
    karray2 = ne.evaluate("kxarray**2 + kyarray**2 + kzarray**2")
    return [kxarray, kyarray, kzarray], karray2