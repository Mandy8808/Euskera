"""SINGLE-FREQUENCY SHOOTING AND EVENT-IDENTIFICATION HELPERS"""
import numpy as np

def shoot(imin: float, imax: float) -> float:
    """
    Compute the midpoint between imin and imax.

    Parameters:
    - imin (float): Lower bound.
    - imax (float): Upper bound.

    Returns:
    - float: The midpoint between imin and imax.
    """
    if imin > imax:
        raise ValueError("imin should be less than or equal to imax")

    return (imin + imax)/2

def freq_shoot(events, nodo, u0, iInterv, rTemp, inv=False):
    """
    Adjusts U bounds based on event crossings.

    Parameters:
    - events: list of arrays, where events[0] and events[1] represent different event crossings.
    - nodo: int, number of expected crossings.
    - u0: float, current u guess.
    - iInterv: tuple (imin, imax), the interval of frequency adjustment.
    - rTemp: float, last recorded event.
    - inv: bool, flag to determine inversion logic.

    Returns:
    updated the u interval and last event crossing.
    -> [umin, umax], rTemp
    """
    imin, imax = iInterv
    events0, events1 = events  # Unpacking for clarity
    i0 = u0  # Extract first value if it's in a list/array

    if events0.size == nodo and events1.size == nodo + 1:
        return [imin, imax], rTemp
    else:
        # Determine whether to increase or decrease the frequency bounds
        target_event = events0 if events0.size > nodo else events1

    if inv:
        imin, imax = (i0, imax) if target_event is events0 else (imin, i0)
    else:
        imin, imax = (imin, i0) if target_event is events0 else (i0, imax)

    rTemp = target_event[-1]  # Last crossing event update

    return [imin, imax], rTemp

def identify(events, nodos, p0Data):
    """
    Identifies which profile's U bound should be modified based on event crossings.

    Parameters:
    - events: list of NumPy arrays, where each index corresponds to a detected event.
    - nodos: list of expected number of crossings for each component.
    - p0Data: list of bools indicating which components are considered in the analysis.

    Returns:
    - sigModif: list of bools indicating which signals should be modified.
    """
    numFields = len(nodos)
    name = [str(i) for i in range(0, 2 * numFields, 2)]
    dicNod = dict(zip(name, nodos))

    # Identify components that are not zeros
    indices = np.array(list(map(int, dicNod.keys())), dtype=int)
    ind = np.array(p0Data, dtype=bool)  # Ensuring a boolean mask
    posit = indices[ind]

    valR = [np.inf] * numFields  # Initialize with infinity
    for i in posit:
        valtemp = []
        nodo = dicNod[str(i)]

        numNod = len(events[i])  # Size of current event
        numdSig = len(events[i+1]) if i+1 < len(events) else 0

        if (numNod == nodo) and (numdSig == nodo + 1):
            valtemp.append(0)
        elif numNod == nodo:
            if numdSig < nodo + 1:
                valtemp.append(events[i+1][nodo-1] if numdSig > 0 else np.inf)
            elif numdSig > nodo + 1:
                valtemp.append(events[i+1][nodo])
        elif numNod > nodo:
            valtemp.append(events[i][nodo])
        else:
            valtemp.append(events[i][-1] if numNod != 0 else 0)

        # Dynamically assign to valR based on index position
        valR[i // 2] = min(valtemp)

    # Determine the signal to modify
    sigModif = [False] * numFields
    test = np.min(valR)
    for i in range(numFields):
        if valR[i] == test:
            sigModif[i] = True
            break

    return sigModif

def freq_shoot2(events, nodo, i0, iInterv, rTemp):
    """
    """
    imin, imax = iInterv[0]
    events = events[0]
    i0 = i0[0]

    if events[0].size == nodo and events[1].size == nodo+1:
        return [imin, imax], rTemp
    elif events[1].size > nodo+1:
        if events[0].size > nodo:  # dos veces por nodo
            imax = i0
            rTemp = events[0][-1]
        else:  # si pasa por cero más veces que 2*nodos se aumenta la w, sino se disminuye
            imin = i0
            rTemp = events[1][-1]
    elif events[1].size <= nodo+1:
        if events[0].size > nodo:  # dos veces por nodo
            imax = i0
            rTemp = events[0][-1]
        else:
            imin = i0
            rTemp = events[1][-1]
    return [imin, imax], rTemp

def identify2(events, nodos, p0Data, info=False):
    """
    """
    dicNod = {'0': nodos[0], '2': nodos[1], '4': nodos[2]}

    # identificando que componentes no son ceros
    # cuando una componente se tomó como cero y se excluye del análisis
    indices = np.fromiter(map(int, dicNod.keys()), dtype=int)
    ind = list(map(bool, p0Data))
    posit = indices[ind]

    valR = [np.inf, np.inf, np.inf]
    for i in posit:
        #if info:
        #    print(events[i], events[i+1])

        valtemp = []
        nodo = dicNod[str(i)]
        numNod = events[i].size; numdSig = events[i+1].size
        if numNod == nodo and numdSig == nodo+1:
            valtemp.append(0)
        elif numNod == nodo:
            if numdSig < nodo+1:
                valtemp.append(events[i+1][nodo-1])
            elif numdSig > nodo+1:
                valtemp.append(events[i+1][nodo])
        elif numNod > nodo:
            valtemp.append(events[i][nodo])
        else:
            if numNod != 0:
                valtemp.append(events[i][-1])
            else:
                valtemp.append(0)

        if i==0:
            valR[0] = min(valtemp)
        elif i==2:
            valR[1] = min(valtemp)
        elif i==4:
            valR[2] = min(valtemp)

    # print(valR)
    sigModif = [False, False, False]
    test = np.min(valR)
    for i in range(3):
        if valR[i]==test:
            sigModif[i]=True
            break

    return sigModif

__all__ = ['shoot', 'freq_shoot', 'identify', 'freq_shoot2', 'identify2']
