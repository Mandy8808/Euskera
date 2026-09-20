"""FREQUENCY ANALYSIS FOR EUSKERA SIMULATIONS"""

import numpy as np

########################################################################################################################
def main_frequency(psi, t, met=2, info=False, tol=1e-10, guess=1.0, save=False):
    """
    Extrae componentes principales de psi y calcula su frecuencia usando uno de dos métodos.

    Parameters:
        psi : np.ndarray
            Arreglo de evolución temporal de estados (1D, 2D, o 3D por ejemplo).
        t : array_like
            Vector de tiempos.
        met : int
            Método a usar (1 o 2).
        info : bool
            Si se desea imprimir información de depuración.
        tol : float
            Tolerancia para considerar una componente como nula.

    Returns:
        data_component : output de frequMet1 o frequMet2
    """
    #ndim = psi.shape[-1] if psi.ndim == 2 else 1  # assuming psi.shape = (n, d) o (n,)
    ndim = len(psi[0])

    psit = [[] for _ in range(ndim)]
    if ndim == 1:  # vector scalar
        for comp in psi:
            psit[0].append(comp[0])
    elif ndim == 2 or ndim == 3:  # vector with d components
        for comp in psi:
            for i in range(ndim):
                psit[i].append(comp[i])
    else:
        raise ValueError(f"psi debe ser un arreglo 1D o 2D/3D, no {ndim}D")

    # taking only the non-null row (vector componets)
    psit = [row for row in psit if not np.all(np.isclose(row, 0, atol=tol))]  # taking only the non-null row (vector componets)

    # calling the choosed method
    if met == 1:
        data_component = frequMet1(t, psit, info=info, guess=guess)
    else:
        data_component = frequMet2(t, psit, info=info, save=save)

    return data_component


########### Frequency
################################################################################

def frequMet1(t_data, psit, guess=1.0, info=False):
    """
    Estimate frequency ω from data assuming:
        ψ(t) = A * exp(ω * t)
    where A is the amplitude and ω is the angular frequency.


    Parameters:
    - t_data: array-like, time values
    - psit: array-like or list of array-like, signal values
    - info: bool, whether to print the regression output

    Returns:
    -
    """
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    # Ensure psi is iterable (list of signals), even if it's just one
    if isinstance(psit[0], (int, float, complex)):
        psit = [psit]

    # Model
    def func(t, A, E):
        return A * np.exp(1j * E * t)

    def funcBoth(t, A, E):
        temp = func(t, A, E)
        return temp.real

    data_component = []
    for comp in psit:
        comp = np.array(comp)
        E_guess = A_guess = guess
        popt, pcov = curve_fit(funcBoth, t_data, comp.real, p0=[A_guess, E_guess])

        if info:
            print(rf"ω estimado: {popt[1]}")
            print(rf"A estimado: {popt[0]}")
            print(rf"Error: ", np.sqrt(np.diag(pcov)))

            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), sharex=False, sharey=False,
                       gridspec_kw=dict(hspace=0.0, wspace=.13))

            ax[0].plot(t_data, comp.real, c='k', label='Real')
            ax[0].plot(t_data, func(t_data, popt[0], popt[1]).real, ls='--', c='r')

            ax[1].plot(t_data, comp.imag, c='k', label='Imaginary')
            ax[1].plot(t_data, func(t_data, popt[0], popt[1]).imag, ls='--', c='r')

            # plt.show()
        data_component.append([popt, pcov])
    return data_component

def frequMet2(t, psit, info=False, save=False):
    """
    Estimate dominant angular frequency ω from signal(s) using FFT.
    Assumes ψ(t) = A * exp(iωt) or similar periodic form.

    Parameters:
    - t: array-like, time vector (assumed uniform spacing)
    - psit: array-like or list of array-like, signal(s)
    - info: bool, print estimated ω for each component
    - save : Boolean or string indicating whether to save the plot

    Returns:
    - List of dominant ω values (rad/s) for each component
    """
    # Compute time step
    dt = t[1] - t[0]

    # Ensure psi is iterable
    if isinstance(psit[0], (int, float, complex)):
        psit = [psit]

    if info:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(nrows=1, ncols=len(psit), figsize=(len(psit) * 5, 4), sharex=False, sharey=False,
                       gridspec_kw=dict(hspace=0.0, wspace=.23))
        if len(psit) == 1:
            ax = [ax]

    data_component = []
    for i, comp in enumerate(psit):
        comp = np.array(comp)

        # FFT (with zero-padding for better resolution)
        nfft = 4 * len(comp)
        psi_fft = np.fft.fft(comp, n=nfft)
        freq = np.fft.fftfreq(nfft, d=dt)  # frequency in Hz (cycles/sec)
        omega = 2 * np.pi * freq  # Convert to angular frequency (rad/s)

        # Magnitude spectrum (normalized)
        magnitude = np.abs(psi_fft) / len(comp)

        if np.iscomplexobj(comp):
            # For complex signals, consider both positive and negative frequencies
            dominant_idx = np.argmax(magnitude)
            omega_dominante = omega[dominant_idx]
            estimated_phase = np.angle(psi_fft[dominant_idx])  # Obtener la fase
        else:
            # For real signals, consider only positive frequencies
            half = len(omega) // 2  # solo para señales reales
            dominant_idx = np.argmax(magnitude[:half])
            omega_dominante = omega[:half][dominant_idx]
            estimated_phase = np.angle(psi_fft[:half][dominant_idx])  # Obtener la fase

        data_component.append(omega_dominante)

        if info:
            print(rf"Frequency ω ≈ {omega_dominante:12.10f} rad/s")
            print(rf"Fase estimada: {estimated_phase:.2f} rad")  # {np.degrees(estimated_phase):.2f} grados


            ax[i].plot(omega, magnitude)
            ax[i].plot(omega, magnitude, marker='o', markersize=2, linestyle='None', c='k', mfc='white')
            ax[i].vlines(omega_dominante, 0, np.max(magnitude), color='r',
                      linestyle='--', lw=1, label=r'$\omega=%5.4f$ rad/s'%omega_dominante)

            ax[i].set_ylim(0, np.max(magnitude)+0.01)
            ax[i].set_xlim(np.min(omega), np.max(omega))
            ax[i].set_title('FFT Magnitude Spectrum', fontsize=12)
            ax[i].set_ylabel('Magnitude')
            ax[i].set_xlabel('Frequency (rad/s)')
            ax[i].legend(loc='upper left', frameon=False, fontsize=12)

    # Save output
    if info:
        if save:
            filename = save if isinstance(save, str) else 'fft_output.pdf'
            fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.tight_layout()
        plt.show()

    return data_component
