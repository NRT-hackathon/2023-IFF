# Example Diffusion Simulation Scripts

These are example equilibration and production scripts used to run a diffusion simulation for codeine. Automated generation of these scripts was not completed, so the following changes need to be manually made in order to use these scripts for different molecules.

### In equil.lmpin:
+ The **create_box** command needs to be modified to include the required number of atom types, bonds, and angles. Generally, you will need 3 + N atom types, 1 + B bond types, and 1 + A angle types, where N, B, and A are the numbers of atom, bond, and angle types in the fragrance molecule respectively. 
+ The **read_data** command needs to be modified. Replace codeine_data.dat with {chemical_name}_data.dat. If you are changing the length of the polyurea chains, PU_20.dat will need to be changed to PU_{number_of_beads_per_chain}.dat
+ If you are changing the simulation box size, **fix 2** needs to be modified accordingly. 
+ To generate additional trials, we varied the seed used in the **velocity** command.

### In production.lmpin:
+ The **group** command needs to be modified to include all atom types that make up the fragrance molecule. Generally, this is 4 to (3+N), where N is the number of atom types in the fragrance molecule. 
