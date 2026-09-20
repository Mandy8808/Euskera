"""POTENTIAL FUNCTIONS FOR SCHRÖDINGER-POISSON SYSTEM"""

import numpy as np
import numexpr as ne

try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftw = None
    pyfftwOpt = False


#############################################################################

########### Potential Functions
#############################################################################
def Upotential(field_components, rho_i, distarray, rkarray2, num_threads, cmass=0, resol=128, obj=None, check=True):
    """
    Compute the initial value of the potential in the Schrödinger-Poisson system.

    This function calculates the potential `φ(x, t)` using the Poisson equation:

        ∇²φ = 4πG m (|ψ|² - <|ψ|²>)

    The potential field is computed via inverse Fourier transform:

        φ(x, t + h) = F⁻¹ [ (-4π / k²) * F(Sum_j |ψ_j(x, tᵢ)|²) ]

    Parameters:
        rho (ndarray): Density field (|ψ|²).
        distarray (ndarray): Distance array for potential adjustment.
        rkarray2 (ndarray): k² = kx² + ky² + kz²
        num_threads (int): Number of threads for parallel computation.
        cmass (float, optional): Central mass contribution (default 0).

    Notice that:
    Sum_j |ψ_j(x, tᵢ)|² = ψ_0(x, tᵢ) + ψ_1(x, tᵢ) + ..

    Returns:
        tuple: (Potential field, k-space grid components)
    """

    ######## Validate inputs
    if check:
        if not isinstance(rho_i, np.ndarray) or not isinstance(distarray, np.ndarray) or not isinstance(rkarray2, np.ndarray):
            raise TypeError("The elements `rho`, `distarray` or  `rkarray2 ` must be NumPy arrays.")

        if rho_i.shape != (field_components, resol, resol, resol):
            raise ValueError(f"`rho` must have shape ({field_components}, {resol}, {resol}, {resol}).")
    ################################################################################################################################

    use_cached_fftw = obj is not None and pyfftwOpt
    if use_cached_fftw:
        rfft_rho, irfft_phi = obj

    # Compute FFT of rho = Sum_j |ψ_j(x, tᵢ)|²
    rho = np.sum(rho_i, axis=0)  # ne.evaluate("sum(rho_i, axis=0)")
    if not use_cached_fftw:
        if pyfftwOpt:
            rfft_rho = pyfftw.builders.rfftn(rho, axes=(0, 1, 2), threads=num_threads)   # Return a pyfftw.FFTW object representing an n-D real FFT
            phik = rfft_rho(rho)  # Compute the N-dimensional discrete Fourier Transform for real input rho=|psi|^2 (i.e. Fourier(rho))
        else:
            phik = np.fft.rfftn(rho, axes=(0, 1, 2))
    else:
        phik = rfft_rho(rho)

    # Compute potential in k-space: φ_k = (-4π / k²) * ρ_k where ρ_k=Fourier(rho)
    with np.errstate(divide='ignore', invalid='ignore'):
        phik = ne.evaluate("-4 * pi * phik / rkarray2", local_dict={'rkarray2': rkarray2,
                                     'phik': phik, 'pi': np.pi})

    phik[0, 0, 0] = 0  # Set the k=0 Fourier mode to zero prior to the final inverse Fourier transform
                       # This imply there is no need to subtract the global average density <|\psi|^2>.

    # Inverse FFT to obtain real-space potential
    if not use_cached_fftw:
        # Initialize potential field
        if pyfftwOpt:
            phisp = pyfftw.zeros_aligned((field_components, resol, resol, resol), dtype='float64')
            irfft_phi = pyfftw.builders.irfftn(phik, axes=(0, 1, 2), threads=num_threads)   # Return a pyfftw.FFTW object representing an n-D real inverse FFT.
            phisp = irfft_phi(phik)   # Compute the N-dimensional discrete inverse FFT for real inputphik (i.e.  F^{-1} (-/k^2) F 4pi |psi(\vec{x}, t_i)|^2)
        else:
            phisp = np.fft.irfftn(phik, s=rho.shape, axes=(0, 1, 2))
    else:
        phisp = irfft_phi(phik)

    # Adjust potential using central mass, avoiding division by zero
    phisp = ne.evaluate("phisp - (cmass / distarray)",
                        global_dict={'distarray': np.where(distarray == 0, np.inf, distarray),
                                     'cmass': cmass, 'phisp': phisp})

    fft_objects = (rfft_rho, irfft_phi) if pyfftwOpt else (None, None)
    return (phisp, rho) if obj is not None else (phisp, rho, fft_objects)