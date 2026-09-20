"""Spectral eigensolver implementation.
Details in https://arxiv.org/pdf/2512.04376
"""

import numpy as np
from scipy.linalg import eig
from .operators import backgroundOper
from .blocks import linBlock, circBlock, radBlock, multBlock
from euskera.tools.tools import progressbar

def espectro(datFunc, util, lamV, alp, polariz, info=False, fplot=False):
   r"""Computes the spectrum of a system.

   Type of polarizations:
   i. $\gamma=0$, $\alpha=0$ if the polarization is linear,
   ii. $\gamma=0$, $\alpha=1$ if the polarization is circular, and
   iii. $\gamma=1$, $\alpha=0$ if the polarization is radial.
   """
   ############################################################################
   # Cheking the in put parameters:
   # Check if the polarization is valid
   if polariz not in ['linear', 'circular', 'radial', 'multifrequency']:
      raise ValueError("Choose a valid polarization: radial, circular, linear, or multifrequency")
   # Check if the parameters are valid
   if polariz == 'linear' and alp != 0:
      raise ValueError("Alpha should be zero for linear polarization.")
   if polariz == 'circular' and alp == 0:
      raise ValueError("Alpha should be non-zero for circular polarization.")
   if polariz == 'radial' and alp != 0:
      raise ValueError("Alpha should be zero for radial polarization.")
   if polariz == 'multifrequency' and alp != 0:
      raise ValueError("Alpha should be zero for multifrequency polarization.")
   # Check if the data function is valid
   if np.any([not callable(i[0]) for i in datFunc]):
      raise ValueError("datFunc should be a list of callable functions.")
   # Check if the utility parameters are valid
   if not isinstance(util, (list, tuple)) or len(util) != 3:
      raise ValueError("util should be a list or tuple of length 3.")
   if not isinstance(lamV, (list, tuple)) or len(lamV) != 2:
      raise ValueError("lamV should be a list or tuple of length 2.")
   if not isinstance(alp, (int, float)):
      raise ValueError("alp should be an integer or float.")
   if not isinstance(info, bool):
      raise ValueError("info should be a boolean value.")
   if not isinstance(fplot, bool):
      raise ValueError("fplot should be a boolean value.")
   ############################################################################

   # Global variables
   J, Nptos, _ = util

   # Background operators
   Sigma0, Ueff, TrianJInv, Rmatriz, D2i_chev, r_dis2, lamStar = backgroundOper(datFunc, util, lamV, alp, info=info)

   # Define the (M_J)_ij matrix
   num = 3  # Block number: Dim(OM_chev) = num*(Nptos-1) x num*(Nptos-1)
   OM_chev = np.zeros((num*2*(Nptos-1), num*2*(Nptos-1)), dtype=complex)

   # Block matrices functions
   blockPola = {
        'linear': linBlock,
        'circular': circBlock,
        'radial': radBlock,
        'multifrequency': multBlock
   }

   # Index patterns for different polarizations
   index_patterns = {
      "linear": zip([0, 0, 1, 1, 2], [0, 1, 0, 1, 2]),
      "circular": zip([0, 0, 1, 1, 2], [0, 1, 0, 1, 2]),
      "radial": zip([0, 0, 1, 1, 2], [0, 1, 0, 1, 2]),
      "multifrequency": zip([0, 0, 0, 1, 1, 1, 2, 2, 2],
                            [0, 1, 2, 0, 1, 2, 0, 1, 2],),
      }

   # Validate polarization
   if polariz not in blockPola:
      raise ValueError("Choose a valid polarization: radial, circular, linear, or multifrequency")

   # Compute block data
   extras = None if polariz == "linear" else alp if polariz == "circular" else [Rmatriz, J]
   data = blockPola[polariz](Nptos, Sigma0, Ueff, TrianJInv, D2i_chev, lamV, lamStar, extras)

   # Fill OM_chev matrix
   for ind, (row, col) in enumerate(index_patterns[polariz]):
      OM_chev[2 * row * (Nptos - 1): 2 * (row + 1) * (Nptos - 1),
              2 * col * (Nptos - 1): 2 * (col + 1) * (Nptos - 1)] = data[ind]

   # Plot if required
   if fplot:
      from euskera.visualization.spectral_plot import plotImag
      for Sigma0i, Ueffi in zip(Sigma0, Ueff):
         plotImag(r_dis2, Sigma0i, Ueffi)  # + J * (J + 1) * Rmatriz

   # Compute eigenvalues and eigenvectors using scipy.linalg.eig
   # Mij@X = lambda X
   # print("===> ", np.sum(OM_chev), "size ", OM_chev.shape)
   eig_values, eig_vectors = eig(OM_chev)

   # Verify eigenvalues and eigenvectors if info is True
   # Valify that the system Ax=Lx
   if info:
      test = [np.allclose(OM_chev @ eig_vectors[:, i] - (eig_values[i] * eig_vectors[:, i]),
                          np.zeros((2 * num * (Nptos - 1)), dtype=complex))
              for i in range(len(eig_values))]
      print("Verifying that Ax = λx holds ->\n", np.array(test))

   # Compute the real Lambda values
   Lambda = 1j * eig_values  # Lambda = -i LambdaReal -> LambdaReal = i Lambda

   # Sort eigenvalues from smallest to largest
   sorted_indices = np.argsort(np.abs(Lambda))
   LambdaSorted = Lambda[sorted_indices]
   eig_vectors_sorted = eig_vectors[:, sorted_indices]

   return LambdaSorted, Lambda, eig_vectors_sorted, r_dis2

def LamJval(datFunc, rMax, lamV, alp, polariz,
            Nptos=500, Jval=[0, 1, 2, 3], real=True,
            info=False, fplot=False):
  """
  Computing the spectro for a set of J-values
  """
  if info:
    print('Computing the spectro for the polarizacion ', polariz)

  # Spectro
  lambdaVal, dataExtra = [], []
  for ind, J in enumerate(Jval):
     util = [J, Nptos, rMax]
     LambdaSorted, _, eig_vectors_sorted, r_dis2 = espectro(datFunc, util, lamV, alp, polariz, info=info, fplot=fplot)

     if real:
         jj = np.real(LambdaSorted) != 0  # reales
         lambReal = LambdaSorted[jj]
         vectReal = eig_vectors_sorted[:, jj]
         lambdaVal.append([J, lambReal])
         dataExtra.append([J, [r_dis2, lambReal, vectReal]])
     else:
         lambdaVal.append([J, LambdaSorted])
         dataExtra.append([J, [r_dis2, LambdaSorted, eig_vectors_sorted]])

     progressbar(ind, len(Jval)-1)
  return lambdaVal, dataExtra
