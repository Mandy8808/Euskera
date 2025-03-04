# scrip solitons 

import numpy as np
import multiprocessing
import euskera as eu
import matplotlib.pyplot as plt

#############################################################
### Referential Profiles
# n = 0
En0 = 2.4538867992147075  # energy eingevalue
Mn0 = 5.491059312232425/np.sqrt(2)
temp = np.load('./Soliton Profile Files/initial_f_n0.npz')
fn0 = temp['arr_0']

# n = 1
En1 = 2.2967598868951598
Mn1 = 12.212985738200105/np.sqrt(2)
temp = np.load('./Soliton Profile Files/initial_f_n1.npz')
fn1 = temp['arr_0']


#############################################################
#### identifying the cpu_processing number
num_threads = multiprocessing.cpu_count()  


#### Symulation Parameters
simulation_parameters = {
"lambda_value": 0,  # Self-interaction coupling constant.
"num_threads": num_threads,  # Set the number of threads to use.
"gridlength": 20,  # Length of the grid in codec units.
"resol": 250,  # Grid resolution. The grid is computed using:
               # linspace(-gridlength/2.0 + gridlength/(2.0*resol), 
               #           gridlength/2.0 - gridlength/(2.0*resol), resol).
"step_factor": 1.,  # Increase this if velocities are low enough that the timestep constraint can be relaxed.
"t0": 0,  # Initial time.
"tmax": 100,  # Total simulation time.
"Plim": 5.6,  # Boundary radius of the soliton.
"rmax": 7.6,  # Two rmax define the minimum separation between solitons 
              # to ensure they are not overlapping.
"cmass": 0,  # Central mass, given in the same units as the soliton mass.
"plott0": False,  # Show the initial profile if True.
"Boverlap": False,  # Check for possible overlap if True.
"info": False  # Print extra information if True.
}

#### Save data information/parameters
salva_data = {
  "format": "npz",  # Format used to save data. Supported formats: "npz" and "hdf5".
  "address": "Data_cuadrado",  # Directory where the data will be saved.
  "save_number": 80,  # Number of frames to save within the time interval [t0, tmax].
  "data_save" : {"grid": True,  # Save grid data.
                 "save_rho": False,   # Save density field (ρ) if True.
                 "save_psi": False,  # Save wavefunction (ψ) if True.
                 "save_phi": True,  # Save potential (φ) if True.
                 "save_plane": True,  # Save plane projections  [:, :, resol/2,]] if True.
                 "save_energies": None,  # Save energy values if specified.
                 "save_line": True}  # Save 1D profile along a specific line [:, resol/2,, resol/2,]]  if True.
}


#############################################################
# Example: Two scalar solitons 

###### Used Model
model = "soliton"

###### Parameters for Components
compfield = 1  # Number of soliton components (for a scalar field: 1, for a vector field: 2, e.g., Ψ = [ψ_x, ψ_y], etc.)

# Component Profiles
sol1 = [fn0]  # If Ψ is a vector (e.g., Ψ = [ψ_x, ψ_y]), specify [fn0, fn0].
sol2 = [fn0]
sol3 = [fn0]
sol4 = [fn0]

# Energy Eigenvalues for Each Component
beta1 = [En0]  # If Ψ is a vector, specify [En0, En0].
beta2 = [En0]
beta3 = [En0]
beta4 = [En0]

###### Configuration Parameters
# Central Positions (x0, y0, z0)
csol1 = [-4., 4., 0.]
csol2 = [-4., -4., 0.]
csol3 = [4., 4., 0.]
csol4 = [4., -4., 0.]

# Velocities
vsol1 = [0., 0., 0.]
vsol2 = [0., 0., 0.]
vsol3 = [0., 0., 0.]
vsol4 = [0., 0., 0.]

# Phases
phase1 = [0] 
phase2 = [0]
phase3 = [0]
phase4 = [0]

# Alpha values (only for the non-self-interacting case; otherwise, set to 1)
alpha = lambda massC, mass: (massC / mass) ** 2  # Conversion factor
alpha1 = [alpha(8, Mn0)]
alpha2 = [alpha(8, Mn0)]
alpha3 = [alpha(8, Mn0)]
alpha4 = [alpha(8, Mn0)]

# Step size (dr) used to computed the radial profiles previusly
delta_x = 0.00001
dr1 = [delta_x]
dr2 = [delta_x]
dr3 = [delta_x]
dr4 = [delta_x]

#############################################################

parameters_sol = [
    {
        "profiles": sol1,
        "positions": csol1,
        "velocities": vsol1,
        "betas": beta1,
        "phases": phase1,
        "alphas": alpha1,
        "dr": dr1
    },
    {
        "profiles": sol2,
        "positions": csol2,
        "velocities": vsol2,
        "betas": beta2,
        "phases": phase2,
        "alphas": alpha2,
        "dr": dr2
    },
    {
        "profiles": sol3,
        "positions": csol3,
        "velocities": vsol3,
        "betas": beta3,
        "phases": phase3,
        "alphas": alpha3,
        "dr": dr3
    },
    {
        "profiles": sol4,
        "positions": csol4,
        "velocities": vsol4,
        "betas": beta4,
        "phases": phase4,
        "alphas": alpha4,
        "dr": dr4
    }
]

# Beginning the simulation
eu.evolve(
    model,
    model_parameters=parameters_sol,
    field_components=compfield,
    salva_data_update=salva_data,
    simulation_parameters_update=simulation_parameters,
    info=False
)


# Generating video
address = "Data_cuadrado"
figData = plt.subplots(nrows=1, ncols=2, figsize=(8, 3.5), sharex=False, sharey=False,
                       gridspec_kw=dict(hspace=0.0, wspace=.13))

dvideo = eu.Visualization(address, figData)  # creando el objeto

n0 = 0  # begin frame
nameP = "save_line"
name2 = "save_plane"
struc = ['-', 2, '#f0784d']
nameV = 'Caso_cuadrado_config'
coord = ['x']
cmapint=[(0, '#050505'),
         (0.05, '#f59876'),
         (0.5, '#f57d51'),
         (0.75, '#f76834'),
          (1, '#f54a0c')]
vconf = [4, 1000000, ['-vcodec', 'libx264']]

text =  '' #r'$n=0, \lambda=0$.'

video = dvideo.video(nameP, struc, nameV, name2=name2, cmapint=cmapint,
                     coord=coord, show=True, n0=n0, plot2D=True, plot3D=False,
                     vconf=vconf, save=True,
                     xlim=(-1, 1), ylim=(-1, 1), text=text)