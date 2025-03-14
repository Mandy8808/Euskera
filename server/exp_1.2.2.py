# scrip solitons 
import os
import sys
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

# Get the parent directory dynamically
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)

import euskera as eu

#############################################################
### Referential Profiles
# n = 0
En0 = 2.4538867992147075  # energy eingevalue
Mn0 = 5.491059312232425/np.sqrt(2)
fn0 = np.load('./Soliton Profile Files/initial_f_n0.npz')['arr_0']

# n = 1
En1 = 2.2967598868951598
Mn1 = 12.212985738200105/np.sqrt(2)
fn1 = np.load('./Soliton Profile Files/initial_f_n1.npz')['arr_0']

#############################################################
#### identifying the cpu_processing number
num_threads = 32 # multiprocessing.cpu_count()  


#### Symulation Parameters
simulation_parameters = {
"lambda_value": 0,  # Self-interaction coupling constant.
"num_threads": num_threads,  # Set the number of threads to use.
"gridlength": 60,  # Length of the grid in codec units.
"resol": 250,  # Grid resolution. The grid is computed using:
               # linspace(-gridlength/2.0 + gridlength/(2.0*resol), 
               #           gridlength/2.0 - gridlength/(2.0*resol), resol).
"step_factor": 1.,  # Increase this if velocities are low enough that the timestep constraint can be relaxed.
"t0": 0,  # Initial time.
"tmax": 200,  # Total simulation time.
"Plim": 8.6,  # Boundary radius of the soliton.
"rmax": 10.6,  # Two rmax define the minimum separation between solitons 
              # to ensure they are not overlapping.
"cmass": 0,  # Central mass, given in the same units as the soliton mass.
"plott0": False,  # Show the initial profile if True.
"Boverlap": False,  # Check for possible overlap if True.
"methodEnerg": 1,  # methodology used to computed the kinetic contribution
"info": False  # Print extra information if True.
}

#### Save data information/parameters
salva_data = {
  "format": "npz",  # Format used to save data. Supported formats: "npz" and "hdf5".
  "address": "Data_122",  # Directory where the data will be saved.
  "save_number": 200,  # Number of frames to save within the time interval [t0, tmax].
  "data_save" : {"grid": True,  # Save grid data.
                 "save_rho": False,   # Save density field (ρ) if True.
                 "save_psi": False,  # Save wavefunction (ψ) if True.
                 "save_phi": False,  # Save potential (φ) if True.
                 "save_plane": True,  # Save plane projections  [:, :, resol/2,]] if True.
                 "save_energies": True,  # Save energy values if specified.
                 "save_line": True}  # Save 1D profile along a specific line [:, resol/2,, resol/2,]]  if True.
}


#############################################################
# Example: 1.2.1

###### Used Model
model = "soliton"
nume_soliton = 2

###### Parameters for Components
compfield = 1  # Number of soliton components (for a scalar field: 1, for a vector field: 2, e.g., Ψ = [ψ_x, ψ_y], etc.)

# Component Profiles
sol1 = [fn1]  # If Ψ is a vector (e.g., Ψ = [ψ_x, ψ_y]), specify [fn0, fn0].
sol2 = [fn0]
solT = [sol1, sol2]

# Energy Eigenvalues for Each Component
beta1 = [En1]  # If Ψ is a vector, specify [En0, En0].
beta2 = [En0]
betaT = [beta1, beta2]

###### Configuration Parameters
# Central Positions (x0, y0, z0)
csol1 = [0, 0, 0]
csol2 = [4, 0, 0]
csolT = [csol1, csol2]

# Velocities
vsol1 = [0., 0., 0.]
vsolT = [vsol1]*nume_soliton

# Phases
phase1 = [0.] 
phaseT = [phase1]*nume_soliton

# Alpha values (only for the non-self-interacting case; otherwise, set to 1)
alpha = lambda massC, mass: (massC / mass) ** 2  # Conversion factor
alpha1 = [alpha(8, Mn0)]
alpha2 = [alpha(0.5, Mn0)]
alphaT = [alpha1, alpha2]

# Step size (dr) used to computed the radial profiles previusly
delta_r = 0.00001
dr1 = [delta_r]
drT = [dr1]*nume_soliton

#############################################################

parameters_sol = {
  model: [{"profiles": solT[i], "positions": csolT[i], "velocities": vsolT[i], "betas": betaT[i],
            "phases": phaseT[i], "alphas": alphaT[i], "dr": drT[i]} for i in range(0, nume_soliton)]
}

# Beginning the simulation
eu.evolve(model_parameters=parameters_sol,
           field_components=compfield,
           salva_data_update=salva_data,
           simulation_parameters_update=simulation_parameters,
           info=False)


# Generating video
#address = "Data_n1_n0_lam_0"
#figData = plt.subplots(nrows=1, ncols=2, figsize=(8, 3.5), sharex=False, sharey=False,
#                       gridspec_kw=dict(hspace=0.0, wspace=.13))

#dvideo = eu.Visualization(address, figData)  # creando el objeto

#n0 = 0  # begin frame
#nameP = "save_line"
#name2 = "save_plane"
#struc = ['-', 2, '#f0784d']
#nameV = 'Caso_n1_n0_lam_0'
#coord = ['x']
#cmapint=[(0, '#050505'),
#         (0.05, '#f59876'),
#         (0.5, '#f57d51'),
#         (0.75, '#f76834'),
#          (1, '#f54a0c')]
#vconf = [4, 1000000, ['-vcodec', 'libx264']]

#text =  '' #r'$n=0, \lambda=0$.'

#video = dvideo.video(nameP, struc, nameV, name2=name2, cmapint=cmapint,
#                     coord=coord, show=True, n0=n0, plot2D=True, plot3D=False,
#                     vconf=vconf, save=True,
#                     xlim=(-1, 1), ylim=(-1, 1), text=text)
