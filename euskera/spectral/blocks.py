"""
Spectral Polarizations blocks implementation.
Details of the linear, circular, and radial block constructions appear in
https://arxiv.org/pdf/2512.04376
"""

import numpy as np

def linBlock(Nptos, Sigma0, Ueff, TrianJInv, D2i_chev, lamV, lamStar, extras=None):
   """
   Constructs the matrix blocks M11, M12, M21, M22, and M33 using the given background operators.
   These matrices appear in the linearized equations (Eq. A2).

   Parameters:
   - Nptos: Number of discretization points.
   - Sigma0: Background sigma matrix.
   - Ueff: Effective potential matrix.
   - TrianJInv: Inverted triangular operator.
   - D2i_chev: Second derivative operator matrix.
   - lamV: Tuple (ln, ls) representing the lambda values.

   Returns:
   - M11, M12, M21, M22, M33: Block matrices.
   """
   ln, ls = np.array(lamV) / lamStar
   size = 2 * (Nptos - 1)

   # Initialize complex matrices
   M11 = np.zeros((size, size), dtype=complex)
   M12 = np.zeros((size, size), dtype=complex)
   M21 = np.zeros((size, size), dtype=complex)
   M22 = np.zeros((size, size), dtype=complex)

   Sigma0, Ueff = Sigma0[0], Ueff[0]
   sigSq = Sigma0 @ Sigma0  # Compute Sigma0^2
   sigma_Trian_sigma = Sigma0 @ TrianJInv @ Sigma0  # Precompute term to avoid redundancy

   # Define block positions
   row, col = 0, 1
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - ln * sigSq
   M22[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - (ln + 2 * ls) * sigSq

   row, col = 1, 0
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - 3 * ln * sigSq - 2 * sigma_Trian_sigma
   M22[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - ln * sigSq

   # M33 is identical to M22
   M33 = np.copy(M22)

   return M11, M12, M21, M22, M33

def circBlock(Nptos, Sigma0, Ueff, TrianJInv, D2i_chev, lamV, lamStar, extras):
   """
   Constructs the matrix blocks M11, M12, M21, M22, and M33 using the given background operators.
   These matrices appear in the linearized equations (Eq. A3).

   Parameters:
   - Nptos: Number of discretization points.
   - Sigma0: Background sigma matrix.
   - Ueff: Effective potential matrix.
   - TrianJInv: Inverted triangular operator.
   - D2i_chev: Second derivative operator matrix.
   - lamV: Tuple (ln, ls) representing the lambda values.
   - alph: Scaling factor.
   - lamStar: Scaling parameter for lambda values.

   Returns:
   - M11, M12, M21, M22, M33: Block matrices.
   """
   ln, ls = np.array(lamV) / lamStar

   alph = extras
   size = 2 * (Nptos - 1)
   # Initialize complex matrices
   M11 = np.zeros((size, size), dtype=complex)  # dtype = 'complex_'
   M12 = np.zeros((size, size), dtype=complex)
   M21 = np.zeros((size, size), dtype=complex)
   M22 = np.zeros((size, size), dtype=complex)
   M33 = np.zeros((size, size), dtype=complex)

   Sigma0, Ueff = Sigma0[0], Ueff[0]
   sigSq = Sigma0 @ Sigma0  # Compute Sigma0^2
   sigma_Trian_sigma = Sigma0 @ TrianJInv @ Sigma0  # Precompute term to avoid redundancy

   row, col = 0, 0
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M12[start_r:end_r, start_c:end_c] = 1j * alph * ls * sigSq
   M21[start_r:end_r, start_c:end_c] = -1j * alph * ((ln + 2 * ls) * sigSq + sigma_Trian_sigma)

   row, col = 0, 1
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - ln * sigSq
   M22[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - (2 * ln + ls) * sigSq - sigma_Trian_sigma
   M33[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - (ln + ls) * sigSq

   row, col = 1, 0
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - (2 * ln + ls) * sigSq - sigma_Trian_sigma
   M22[start_r:end_r, start_c:end_c] = D2i_chev + Ueff - ln * sigSq
   M33[start_r:end_r, start_c:end_c] = D2i_chev + Ueff -(ln + ls) * sigSq

   row, col = 1, 1
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M12[start_r:end_r, start_c:end_c] = 1j * alph * ((ln + 2 * ls) * sigSq + sigma_Trian_sigma)
   M21[start_r:end_r, start_c:end_c] = -1j * alph * ls * sigSq

   return M11, M12, M21, M22, M33

def radBlock(Nptos, Sigma0, Ueff, TrianJInv, D2i_chev, lamV, lamStar, extras):
   """
   Constructs the matrix blocks M11, M12, M21, M22, and M33 using the given background operators.
   These matrices appear in the linearized equations (Eq. A4).

   Parameters:
   - Nptos: Number of discretization points.
   - Sigma0: Background sigma matrix.
   - Ueff: Effective potential matrix.
   - TrianJInv: Inverted triangular operator.
   - D2i_chev: Second derivative operator matrix.
   - lamV: Tuple (ln, ls) representing the lambda values.
   - alph: Scaling factor.
   - lamStar: Scaling parameter for lambda values.

   Returns:
   - M11, M12, M21, M22, M33: Block matrices.
   """
   ln, ls = np.array(lamV) / lamStar
   Rmatriz, J = extras

   size = 2 * (Nptos - 1)
   # Initialize complex matrices
   M11 = np.zeros((size, size), dtype=complex)  # dtype = 'complex_'
   M12 = np.zeros((size, size), dtype=complex)
   M22 = np.zeros((size, size), dtype=complex)

   Sigma0, Ueff = Sigma0[0], Ueff[0]
   sigSq = Sigma0 @ Sigma0
   sigma_Trian_sigma = Sigma0 @ TrianJInv @ Sigma0  # Precompute term to avoid redundancy

   row, col = 0, 1
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r: end_r, start_c: end_c] = D2i_chev + Ueff - 2 * Rmatriz - ln*sigSq
   M12[start_r: end_r, start_c: end_c] = 2 * np.sqrt(J * (J + 1)) * Rmatriz
   M22[start_r: end_r, start_c: end_c] = D2i_chev + Ueff - (ln + 2 * ls) * sigSq

   row, col = 1, 0
   start_r, end_r = row * (Nptos - 1), (row + 1) * (Nptos - 1)
   start_c, end_c = col * (Nptos - 1), (col + 1) * (Nptos - 1)
   M11[start_r: end_r, start_c: end_c] = D2i_chev + Ueff - 2 * Rmatriz - 3 * ln * sigSq - 2 * sigma_Trian_sigma
   M12[start_r: end_r, start_c: end_c] = 2 * np.sqrt(J * (J + 1)) * Rmatriz
   M22[start_r: end_r, start_c: end_c] = D2i_chev + Ueff - ln * sigSq

   M21 = np.copy(M12)
   M33 = np.copy(M22)

   return M11, M12, M21, M22, M33

def MultMii(Nptos, Sigma0i, TSigma0_Sq, Ueffi, TrianJInv, D2i_chev, ln):
   """
   Mii -> multifrequency
   """
   size = 2 * (Nptos - 1)
   Mii = np.zeros((size, size), dtype=complex)  # Initialize complex matrices

   sigSq = Sigma0i@Sigma0i
   sigma_Trian_sigma = Sigma0i @ TrianJInv @ Sigma0i

   # Block (0, 1)
   Mii[0*(Nptos - 1): 1*(Nptos - 1), 1*(Nptos - 1): 2*(Nptos - 1)] = D2i_chev + Ueffi - ln * TSigma0_Sq

   # Block (1,0)
   Mii[1*(Nptos - 1): 2*(Nptos - 1), 0*(Nptos - 1): 1*(Nptos - 1)] = (
      D2i_chev + Ueffi - ln * (2 * sigSq + TSigma0_Sq) - 2 * sigma_Trian_sigma
   )

   return Mii

def MultMij(Nptos, Sigma0i, Sigma0j, TrianJInv, ln):
   """
   Mij -> multifrequency
   """
   size = 2 * (Nptos - 1)
   Mij = np.zeros((size, size), dtype=complex)  # Initialize complex matrices

   sigSqij = Sigma0i@Sigma0j
   sigmai_Trian_sigmaj = Sigma0i @ TrianJInv @ Sigma0j

   # Block (1, 0)
   Mij[1*(Nptos - 1): 2*(Nptos - 1), 0*(Nptos - 1): 1*(Nptos - 1)] = -2 * (ln * sigSqij + sigmai_Trian_sigmaj)

   return Mij

def multBlock(Nptos, Sigma0, Ueff, TrianJInv, D2i_chev, lamV, lamStar, extras=None):
   """
   Construct a block matrix using `MultMii` and `MultMij` functions.
   """
   ln, _ = np.array(lamV) / lamStar

   TSigma0_Sq = sum(Sigma0i @ Sigma0i for Sigma0i in Sigma0)  # summation

   Mij = []
   for i in range(3):
      Sigma0i, Ueffi = Sigma0[i], Ueff[i]
      for j in range(3):
         Sigma0j = Sigma0[j]
         if j != i:
            Mij.append(MultMij(Nptos, Sigma0i, Sigma0j, TrianJInv, ln))
         else:
            Mij.append(MultMii(Nptos, Sigma0i, TSigma0_Sq, Ueffi, TrianJInv, D2i_chev, ln))

   return Mij  # M11, M12, M13, M21, M22, M23, M31, M32, M33
