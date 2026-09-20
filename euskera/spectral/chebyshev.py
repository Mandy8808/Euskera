"""Spectral chebyshev implementation.
Details in https://arxiv.org/pdf/2512.04376
"""

import numpy as np

def cheb(op):
    '''Chebyshev polynomial differentiation matrix.
    Ref.: Trefethen's 'Spectral Methods in MATLAB' book.
    N - size of diff matrix - op+1 where op is polynomial order.
    '''
    N = op + 1
    if N == 1:
        return np.array([[0.0]]), np.array([1.0])  # Edge case for low order

    # Creating the Chebyshev points
    j = np.arange(N)
    x = np.cos(j * np.pi / (N - 1))

    # Weights for differentiation matrix
    c = np.ones(N)
    c[0], c[-1] = 2.0, 2.0  # Special weights for first and last points
    c *= (-1.0)**j
    c = c.reshape(N, 1)

    # Compute matrix entries
    X = np.tile(x.reshape(N, 1), (1, N))
    dX = X - X.T  # Difference matrix

    D = np.dot(c, (1.0/c).T) / (dX + np.eye(N))  # Compute off-diagonal elements
    np.fill_diagonal(D, 0)  # Zero diagonal elements first

    # Compute diagonal elements
    D[np.diag_indices(N)] = -D.sum(axis=1)

    return D, x

def cheb2(op):
    '''Chebushev polynomial differentiation matrix.
       Ref.: Trefethen's 'Spectral Methods in MATLAB' book.
       N - size of diff matrix - op+1 where op is polynomial order.
    '''

    N = op + 1
    # creating the Chebyshev points
    j = np.arange(N)
    x = np.cos(j*np.pi/(N-1))

    #  off-diagonal entries -> Dij
    Dtemp = np.ones(N)
    Dtemp[0], Dtemp[N-1] = 2.0, 2.0
    Dtemp *= (-1.0)**j
    Dtemp = Dtemp.reshape(N, 1)
    Dtemp = np.dot(Dtemp, (1.0/Dtemp).T)

    X = np.tile(x.reshape(N, 1), (1, N))
    dX = X - X.T  # distance difference (dX is a matrix)

    Dtemp = Dtemp/(dX+np.eye(N))  # eye Return a N+1 array with ones on the diagonal and zeros elsewhere.

    # diagonal entries
    Dii = np.diag(Dtemp.sum(axis=1))  # sum by the row, and construct a diagonal array

    # Dn
    D = Dtemp - Dii

    return D, x

def chevQuant(util, info=False):
   """
   Map between [-1, 1] -> [0, rMax]
   -> radii
   -> differentiation matrix
   """
   # Unpacking tuple
   _, Nptos, rMax = util

   # Compute the Chebyshev differentiation matrix
   # D_chev, x_chev = cheb(Nptos - 1)  # Nptos = order + 1, so pass Nptos-1
   D_chev, x_chev = cheb2(Nptos)  # Nptos = order + 1, so pass Nptos-1

   # Scaling from [-1, 1] to [0, rMax] using x_chev = 2(r/rMax) - 1
   r_dis = (1. - x_chev) * (rMax / 2.)

   if info:
      # Checking rescaling
      print('Checking the distance scaling:', np.isclose(r_dis[0], 0.), np.isclose(r_dis[-1], rMax))

   return r_dis, D_chev, x_chev
