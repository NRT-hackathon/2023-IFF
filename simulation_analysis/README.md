# Simulation Analysis Scripts

## Diffusion Coefficient
The msd.py script is used to compute the diffusion coefficient from the diffusion coefficient simulation trajectories. Specificially, this will load all files ending with '.lammpstrj' in the current directory and compute the diffusion coefficient. Note that this should be used in conjunction with the 'diffusing_molecs.lammpstrj' output from the diffusion simulations (not 'production.lammpstrj', which also includes the polyurea chains for visualization, etc.)

The script plot loglog plots of the mean squared displacement (msd) vs. simulation time,  and print the mean and standard deviation the diffusion coefficient computed across all the trajectories in the directory.

**To use this code for different fragrance molecules** the mass dictionary must be modified to include the bead types included in the molecule and their respective masses (the bead types are the same as those in lammps, so here they start with 4).

## Partition Coefficient
Well, actually we only compute the free energy change with the TI_integral.py script, and the partition coefficient requires an additional calculation. The TI_integral.py script integrates the 'fdti10.fep' files that are output from the partition simulations. This will loop through all files ending with '.fep' in the current directory and, plot the derivative data, and compute the change in free energy using the trapezoid method of integration over the computed derivatives. This will compute the mean and standard deviation of free energy changes in all files in the directory.

The free energy change we compute here is the free energy change associated with removing the fragrance molcule from its environment (water, polyurea, or more fragrance molecules).
When combining multiple of these free energies, we can calculate the free energy change between the two environments. For example if $\delta G_{water}$ and $\delta G_{pu}$ represent the free energy changes of removing the fragrance molecule from water and polyurea respectively, we can compute the partition coefficient for the fragrance molecule between water and polyurea ($K_{water/pu}$) using:
$log(K_{water/pu}) = \dfrac{\delta G_{water} - \delta G_{pu}}{2.303RT}$ where R is the molar Boltzmann constant (8.314e-3 kJ/mol.K) and T is the temperature.

