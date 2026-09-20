# euskera v1.0
# conserved quantities file
import numpy as np
import numexpr as ne

# Check if pyFFTW is available
try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftw = None
    pyfftwOpt = False

########################################################################################################################
def Conserv(data, comp_conserv, simulation_parameters,
            obj2=None, methodEnerg=1, pyfftwOpt=False,
            max_pos=None):
    """
    Compute conserved quantities in a numerical simulation.
    """
    psi, rho, phisp, distarray, karray2, kvec = data
    kvec2 = [kvec[0].flatten(), kvec[1].flatten(), kvec[2].flatten()]

    # Parameters
    resol = simulation_parameters["resol"]
    gridlength = simulation_parameters["gridlength"]
    num_threads = simulation_parameters["num_threads"]
    ne.set_num_threads(num_threads)

    # FFT object initialization
    if obj2 is None:
        if pyfftwOpt:
            fft_psi = pyfftw.builders.fftn(psi, axes=(1, 2, 3), threads=num_threads)
            ifft_funct = pyfftw.builders.ifftn(psi, axes=(1, 2, 3), threads=num_threads)
        else:
            fft_psi = lambda psi: np.fft.fftn(psi, axes=(1, 2, 3))
            ifft_funct = lambda psi: np.fft.ifftn(psi, axes=(1, 2, 3))

        obj2 = [fft_psi, ifft_funct]

    Vcell = (gridlength / float(resol))**3  # volumen element

    # Define available functions
    comp = {
        "Numb_Part": "Npar(psi, Vcell)", #"Numb_Part": "Npar(rho, Vcell)",
        "Energ": "Energ(rho, psi, phisp, Vcell, distarray, karray2, kvec2, obj2, simulation_parameters, method=methodEnerg)",
        "Pi": "Pi(psi, kvec2, Vcell, obj2)", "Ji": None,
        "Frequency": "psiData(psi, max_pos)"
    }

    indx = indy = indz = None
    # Restrict eval() to known functions
    safe_globals = {"Npar": Npar, "Energ": Energ, "Pi": Pi, "rho": rho, "psi": psi,
                    "phisp": phisp, "Vcell": Vcell, "distarray": distarray,
                    "karray2": karray2, "kvec2": kvec2, "obj2": obj2,
                    "simulation_parameters": simulation_parameters, "methodEnerg": methodEnerg,
                    "psiData": psiData, "max_pos": max_pos
    }


    data_out = [eval(comp[cant], safe_globals) for cant, opt in comp_conserv.items() if opt and cant in comp]

    return np.array(data_out, dtype=object)

########### Particle Number
################################################################################
#def Npar(rho, Vcell):
#    return Vcell * np.sum(rho)

def Npar(psi, Vcell):
    """
    Estimate 'number of particles' or total norm from discretized wavefunction.

    Parameters:
    - psi: list of arrays, each representing a wavefunction component over a 3D grid

    Returns:
    - List of total 'mass' or norm for each component: ∫ |ψ|² dV
    """
    mass_comp = []
    for comp in psi:
        comp = np.array(comp)
        mass = Vcell * np.sum(np.abs(comp)**2)
        mass_comp.append(mass)

    return [sum(mass_comp), mass_comp]

########### Total Energy
################################################################################
def Energ(rho, psi, phisp, Vcell, distarray,
          karray2, kvec, obj, simulation_parameters, method=1):
    """
    Compute the total energy of the system, including:
    - Gravitational potential energy
    - Self-interaction energy (condensate & field)
    - Kinetic energy

    Parameters:
    - rho: Density array
    - psi: Wavefunction array
    - phisp: Scalar potential
    - Vcell: Volume element
    - distarray: Distance array for potential calculations
    - karray2: Squared wavenumber array
    - kvec: Wavenumber components
    - obj: FFT and IFFT objects
    - simulation_parameters: Dictionary containing physical parameters
    - method: Method choice for kinetic energy computation

    Returns:
    - Total energy (float)
    """
    # Compute energy contributions
    contributions = [
        centpotetE(rho, distarray, simulation_parameters),
        selfinterCondensateE(rho, phisp, distarray, simulation_parameters),
        selfinterFieldE(rho, simulation_parameters),
        kintE(psi, karray2, obj, kvec=kvec, method=method)
    ]

    # Compute total energy
    TotEnerg = Vcell * np.sum(contributions)

    return TotEnerg


### Energy Contributions
################################################################################
def centpotetE(rho, distarray, simulation_parameters):
    """
    Compute the gravitational potential energy density.
    """
    # parameters used
    cmass = simulation_parameters["cmass"]
    with np.errstate(divide='ignore', invalid='ignore'):
        dens_centE = ne.evaluate("real((-cmass/distarray) * rho)")

    return np.sum(dens_centE, dtype=np.float64)

def selfinterCondensateE(rho, phisp, distarray, simulation_parameters):
    """
    Compute the self-interaction energy of the condensate.
    """
    # Compute the energy density associated with the self-interaction
    dens_Eup = ne.evaluate("real(0.5 * (phisp * rho))")  # phisp contains both (the central and psi) contributions
    dens_Eup = np.sum(dens_Eup, dtype=np.float64)

    # Gravitational potential energy of the condensate
    dens_centE = centpotetE(rho, distarray, simulation_parameters)

    return dens_Eup - dens_centE

def selfinterFieldE(rho, simulation_parameters):
    """
    Compute the self-interaction energy due to the field.
    """
    # Retrieve lambda_value from the simulation parameters
    lambda_value = simulation_parameters["lambda_value"]

    # Calculate energy density associated with the field self-interaction
    dens_lambd = ne.evaluate("0.25 * lambda_value * rho**2")

    # Return the sum of energy densities
    return np.sum(dens_lambd, dtype=np.float64)

def kintE(psi, karray2, obj, kvec=None, method=1):
    """
    Compute the kinetic energy of the system.
    """
    if method==1:
        fft_psi, ifft_funct = obj

        funct = fft_psi(psi)
        funct = ne.evaluate("karray2 * funct")
        funct = ifft_funct(funct)

        dens_Ekin_i = ne.evaluate("real(0.5 * conj(psi) * funct)")
        dens_Ekin = np.sum(dens_Ekin_i, axis=0)
                    #ne.evaluate("sum(dens_Ekin_i, axis=0)")
    else:
        fft_psi, ifft_funct = obj

        # Compute funct and functconj, store copies to avoid modification
        funct = np.copy(fft_psi(psi))
        psiconj = np.conjugate(psi)
        functconj = np.copy(fft_psi(psiconj))

        # Initialize kinetic energy density array
        dens_Ekin_i = np.zeros_like(psi, dtype='complex128')
        for i in range(3):
            coef1 = np.copy(ifft_funct(kvec[i] * funct))  # Store copy before modifying
            coef2 = np.copy(ifft_funct(kvec[i] * functconj))  # Store copy before modifying
            dens_Ekin_i += coef1 * coef2  # Update the energy density
        dens_Ekin = -0.5 * np.sum(np.real(dens_Ekin_i), axis=0)
                   #-0.5 * ne.evaluate("sum(real(dens_Ekin_i), axis=0)")
    return np.sum(dens_Ekin)
########################################

def Pi(psi, kvec, Vcell, obj):
    fft_psi, ifft_funct = obj
    funct = fft_psi(psi)
    functconj = np.conjugate(funct)
    dens_P_i = []
    for k in range(3):
        temp = ifft_funct(kvec[k] * funct)
        temp = ne.evaluate("sum(functconj * temp, axis=0)")
        dens_P_i.append(temp)

    return [Vcell * np.sum(comp) for comp in dens_P_i]
########################################

def psiData(psi, max_pos):
    data_comp = [comp[ind[0], ind[1], ind[2]] for ind, comp in zip(max_pos, psi)]
    return data_comp