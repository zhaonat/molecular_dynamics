# molecular_dynamics
exploration of MD simulations in MATLAB. A lennard-Jones potential is used for both 2D and 3D (the function is easily differentiable so we can calculate the force and hence acceleration of the particles, which is needed for the Verlet algorithm).

# 2D Simulations
practice for 3D, just want to get the basic idea of the Verlet integrator demonstrated. 

## simulations
contains all the main scripts which can be directly run in the MATLAB window. <br>
For the 2D case, the main verlet algorithms are implemented directly inside the script (instead of partitioned as functions) since they are not too complicated.

![fig](sample_figures/2D_Autocorrelation_Plt_at_T=5.png?raw=true)

 
# 3D simulations

## run_md
This folder contains all the requisite functions required to run the verlet update in the scripts of the simulations folder. Contrary to the 2D case, the extra degrees of freedom in the 3D case makes it a bit more convenient to modularize the Verlet algorithm.

## simulations
Contains all the scripts to execute the molecular dynamics simulation of particles sitting in an Lennard-Jones potential.
![fig](sample_figures/3D_trajectory_sample.png?raw=true)
