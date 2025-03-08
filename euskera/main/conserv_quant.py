# euskera v1.0
# conserved quantities file
import sys
import os
import numpy as np
import numexpr as ne

# Check if pyFFTW is available
try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftw = None
    pyfftwOpt = False


def Conserv(data, cons, field_components, simulation_parameters, obj=None, method=1):
    """ 
    """
    psi, rho, phisp, distarray, karray2, kvec = data
    kvec2 = [kvec[0].flatten(), kvec[1].flatten(), kvec[2].flatten()]
    
    resol = simulation_parameters["resol"]
    gridlength = simulation_parameters["gridlength"]
    num_threads = simulation_parameters["gridlength"]
    
    fft_psi = pyfftw.builders.fftn(psi, axes=(1, 2, 3), threads=num_threads)
    ifft_funct = pyfftw.builders.ifftn(psi, axes=(1, 2, 3), threads=num_threads)
    obj = [fft_psi, ifft_funct] if not obj else obj
    
    Vcell = (gridlength/float(resol))**3
    
    comp = {
        "Numb_Part": Npar(rho, Vcell),
        "Energ": Energ(rho, psi, phisp, Vcell, distarray, field_components,
                       karray2, kvec2, obj, simulation_parameters,
                       method=method),
        "Pi": Pi(psi, kvec2, Vcell, obj), "Ji": None
    }
    
    data_out = []
    for cant, opt in cons.items():
        if opt:
            temp = comp[cant]
            data_out.append(temp)
    return data_out

########################################  
def Npar(rho, Vcell):
    return Vcell * np.sum(rho)

########################################
def Energ(rho, psi, phisp, Vcell, distarray, field_components, 
          karray2, kvec, obj, simulation_parameters, method=1):
    """
    """
    # Gravitational potential energy density associated with the central potential
    contributions = [centpotetE(rho, distarray, simulation_parameters),
                     selfinterCondensateE(rho, phisp, distarray, simulation_parameters),
                     selfinterFieldE(rho, simulation_parameters),
                     kintE(psi, karray2, obj, simulation_parameters, field_components, kvec=kvec, method=method)]
    TotEnerg = Vcell * np.sum(contributions)
    return TotEnerg

def centpotetE(rho, distarray, simulation_parameters):
    cmass, resol = simulation_parameters["cmass"], simulation_parameters["resol"]
    egyarr = pyfftw.zeros_aligned((resol, resol, resol), dtype='float64')
    dens_centE = ne.evaluate("real((-cmass/distarray) * rho)")
    return np.sum(dens_centE)

def selfinterCondensateE(rho, phisp, distarray, simulation_parameters):
    dens_Eup = ne.evaluate("real(0.5 * (phisp * rho))")  # notice that phisp have both,
                                                   # the central and psi contribution
    dens_Eup = np.sum(dens_Eup)
    dens_centE = centpotetE(rho, distarray, simulation_parameters)
    return dens_Eup - dens_centE

def selfinterFieldE(rho, simulation_parameters):
    lambda_value = simulation_parameters["lambda_value"]
    dens_lambd = ne.evaluate("0.25 * lambda_value * rho**2")
    return np.sum(dens_lambd)

def kintE(psi, karray2, obj, simulation_parameters, field_components, kvec=None, method=1):
    if method==1:
        fft_psi, ifft_funct = obj
        funct = fft_psi(psi)
        funct = ne.evaluate("-karray2 * funct")
        funct = ifft_funct(funct)
        dens_Ekin_i = ne.evaluate("real(-0.5 * conj(psi) * funct)")
        dens_Ekin = ne.evaluate("sum(dens_Ekin_i, axis=0)")
    else:
        resol = simulation_parameters["resol"]
        fft_psi, ifft_funct = obj
        funct = fft_psi(psi)
        functconj = np.conjugate(funct)
        dens_Ekin_i = pyfftw.zeros_aligned((field_components, resol, resol, resol), dtype='complex128')
        for i in range(3):
            dens_Ekin_i += ifft_funct(kvec[i]*functconj) * ifft_funct(kvec[i]*funct)
        dens_Ekin = -0.5 * ne.evaluate("sum(dens_Ekin_i, axis=0)")
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