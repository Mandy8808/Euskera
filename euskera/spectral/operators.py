"""Spectral operators implementation."""

import numpy as np
from .chebyshev import chevQuant

def backgroundOper(datFunc, util, lamV, alp, info=False):
   """
   Computes the background operators:
   Sigma0, Ueff, TrianJInv, Rmatriz, D2i_chev, r_dis_inner, lamStar
   """

   # Unpacking the parameters
   J, Nptos, rMax = util

   # Unpacking background quantities
   fsN, fuN = datFunc

   # Computing the \lambda_*
   ln, ls = lamV
   lamStar = abs(ln + alp**2 * ls) if abs(ln + alp**2 * ls) !=0 else 1

   # Mapping [-1, 1] -> [0, rMax]
   utilphy = [J, Nptos, rMax]
   r_dis, D_chev, _ = chevQuant(utilphy, info=info)

   # Creating the matricial operators
   # (NOTICE that the homogeneous Dirichlet boundary conditions are defined by omitting the first and last element)
   r_dis_inner = r_dis[1:Nptos]  # radii values excluding boundaries

   # vector of Matrix Sigma0
   Sigma0 = [np.diag(fs(r_dis_inner)) for fs in fsN]  # [σ_x^0, σ_y^0, ...]

   # temporal matrices
   Rmatriz = np.diag(1 / r_dis_inner**2)  # Rm = 1/r²
   Jmatriz = J * (J + 1) * Rmatriz  # Jm = J(J+1)/r²

   # vector of Matrix Ueff
   Ueff = []
   for fu in fuN:
      U0 = np.diag(fu(r_dis_inner))  # potential u0 = E - Δ⁻¹(|σ₁⁰|²)
      Ueff.append(U0 - Jmatriz)  # Ueff = u0 - J(J+1)/r²

   # Derivative operators
   D2_chev = np.dot(D_chev, D_chev) / ((rMax / 2) ** 2)  # Computing D²

   # Extracting the interior matrix (ignoring boundaries)
   D2i_chev = D2_chev[1:Nptos, 1:Nptos]

   # Operator TrianJInv = [4*D²/r_star² - J(J+1)/r²]⁻¹
   temp = D2i_chev - Jmatriz
   TrianJInv = np.linalg.inv(temp)  # (D² - J(J+1)/r²)⁻¹

   if info:
      print('Dimensions of the matrices Sigma0 and Ueff:', Sigma0[0].shape, Ueff[0].shape)
      print('\nChecking the inverse of D²...')
      check = np.allclose(np.dot(temp, TrianJInv), np.eye(Nptos - 1))  # Check if inverse is correct
      print('Inverse check:', check)

   return Sigma0, Ueff, TrianJInv, Rmatriz, D2i_chev, r_dis_inner, lamStar
