# Release Profile Solver

The solve_release_bootstrap.py file solves the release profile using the macroscale model used by Muro Sune, et al. using the example data from the simulations we ran. The model is sovled for bootstrapped samples of the partition and diffusion coefficients based on normal distributions determined from the mean and standard deviations of the values obtained from simulation.

The solver takes 5 inputs, the parameter list, a list of the standard deviations of each parameter (for bootstrapping), the number of bootstrap samples to use, the maximum time, and the number of time points to solve at. 

To run the release sovler, the parameter list and standard deviation list must be updated following the examples in the script. The number of bootstrap samples (N_bootstraps), the maximum time (t_max), and the number of timepoints (n_times) can all be changed in the bootstrap_solve function as desired.

This can be adapted to run the solver without bootstrapping, simply by setting the standard deviation list to a list of zeros the same length as the parameter array and setting N_bootstraps to 1.

The file currently includes the example data from the two cases that we simulated (all additional parameters, such as the microcapsule sizes, total amount of   from Muro Sune, et al.). 

**Case 1: Codeine**  
+ Diffusion Coefficient: D = 1.954e-9 (8.3e-11) $m^2/s$  
+ Partition Coefficient between polyurea and inside capsule (pure codeine): $log(K_{m/d})$ = 4.25 (0.9)  
+ Partition Coefficient between polyurea and release medium (water): $log(K_{m/d})$ = 12.9 (0.73)  

**Case 2: Hexanal**  
+ Diffusion Coefficient: D = 7.16e-9 (1.86e-9) $m^2/s$  
+ Partition Coefficient between polyurea and inside capsule (pure hexanal): $log(K_{m/d})$ = 2.40 (0.39)  
+ Partition Coefficient between polyurea and release medium (water): $log(K_{m/d})$ = 5.70 (0.38)  

