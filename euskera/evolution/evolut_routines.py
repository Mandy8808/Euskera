"""TIME-EVOLUTION ROUTINES FOR SCHRÖDINGER-POISSON SYSTEM"""

import time
import numpy as np
import numexpr as ne

# Check if pyFFTW is available
try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftw = None
    pyfftwOpt = False

from euskera.observables import simulation_conserv_quant as cq
from euskera.evolution import potential as pt
from euskera.io import save_data as sv
from euskera.tools import tools as to

###################################################################################################

########### Time-evolution function
#############################################################################

def PKP(field_components, fields, param, distDat, num_steps, obj, kvec, simulation_parameters, comp_conserv, data_save_obj, info=False):
    """
    Time evolution of wavefunction psi using the Schrödinger-Poisson system.

    Implements:
        ψ(x, t+h) = exp(-i h φ(x, t+h)/2) * F⁻¹ [ exp(-i h k²/2) * F[exp(-i h φ(x, t)/2) * ψ(x, t)] ]

    Parameters:
        fields (tuple): (phisp, psi), the potential and wavefunction.
        param (tuple): (actual_num_steps, ht, halfstepornot), simulation parameters.
        num_steps (int): Number of evolution steps.
        obj (list): FFTW objects [rfft_rho, irfft_phi].
        pyfftwOpt (bool, optional): Use pyFFTW if available (default: True).

    Returns:
        tuple: Updated (phisp, psi).
    """
    tinit = time.time()

    # Unpack parameters
    num_steps, ht, halfstepornot, its_per_save, num_threads, cmass, resol, lambda_value, max_pos = param
    phisp, psi, rho = fields
    distarray, karray2, rkarray2 = distDat
    methodEnerg = simulation_parameters["methodEnerg"]

    # Select FFT and IFFT functions
    if pyfftwOpt:
        # Using axes=(-3, -2, -1) or  (1, 2, 3) keeps each system (n=field_components) independent because:
	    # 1. It only transforms the internal dimensions of each system (resol, resol, resol).
        # 2. It does not mix the first dimension (field_components), which represents different systems.
        # 3. It is used when each “box” of size (resol, resol, resol) must be analyzed in frequency without influencing the other.
        fft_psi = pyfftw.builders.fftn(psi, axes=(1, 2, 3), threads=num_threads)
        ifft_funct = pyfftw.builders.ifftn(psi, axes=(1, 2, 3), threads=num_threads)
    else:
        fft_psi = lambda psi: np.fft.fftn(psi, axes=(1, 2, 3))
        ifft_funct = lambda psi: np.fft.ifftn(psi, axes=(1, 2, 3))

    # Precompute exponential term for k-space evolution
    exp_k = ne.evaluate("exp(-1j * 0.5 * ht * karray2)")

    # Time evolution loop
    systema0 = "exp(-1j * 0.5 * ht * (0.5 * lambda_value * rho + phisp)) * psi" if lambda_value else "exp(-1j * 0.5 * ht * phisp) * psi"
    systema = "exp(-1j * ht * (0.5 * lambda_value * rho + phisp)) * psi" if lambda_value else "exp(-1j * ht * phisp) * psi"

    ttot, count = 0, 1
    for i in range(int(num_steps)):
        psi, halfstepornot = (ne.evaluate(systema0), False) if halfstepornot else (ne.evaluate(systema), False)

        funct = fft_psi(psi)
        funct = ne.evaluate("funct * exp_k")
        psi = ifft_funct(funct)
        rho_i = ne.evaluate("real(abs(psi)**2)")

        # computing the potential
        (phisp, rho) = pt.Upotential(field_components, rho_i, distarray, rkarray2, num_threads, cmass=cmass, resol=resol, obj=obj,
                              check=False)

        # Save step
        if (((i + 1) % its_per_save) == 0) and halfstepornot == False:
            # Next if statement ensures that an extra half step is performed at each save point
            # psi = ne.evaluate("exp(-1j * 0.5 * ht * phisp) * psi")
            psi = ne.evaluate(systema0)
            rho = ne.evaluate("sum(rho_i, axis=0)")  # rho = np.sum(rho_i, axis=0)
            halfstepornot = True

            # Calculate energies and save data
            data = [psi, rho, phisp, distarray, karray2, kvec]
            cData = cq.Conserv(data, comp_conserv, simulation_parameters, obj2=[fft_psi, ifft_funct], methodEnerg=methodEnerg, max_pos=max_pos)  # method=1
            data = ([None, None, None], rho, psi, phisp, cData)
            sv.fdata_save(ti=count, data=data, data_save_obj=data_save_obj, resol=resol, end=False, physical_time=(i + 1) * ht)
            count += 1

        # Time tracking
        tint = time.time() - tinit
        if info:
            print('cpu time:', tint, '\n')
        ttot += tint
        tinit = time.time()

        # Update progress bar
        to.progressbar(i, num_steps-1, bar_length=20, progress_char='#')

    sv.fdata_save(ti=None, data=([None, None, None], None, None, None, None),
               data_save_obj=data_save_obj, resol=None, end=True)

    print ('\n')
    print("Complete. Total time -> ", ttot)
    print('save_frames', count, '\n')

    return rho, phisp