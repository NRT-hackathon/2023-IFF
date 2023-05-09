# Example Partition Simulation Scripts

These are example equilibration and production scripts used to run partition simulations for hexanal in hexanal (self), polyurea (pu), and water. Automated generation of these scripts was not completed, so the following changes need to be manually made in order to use these scripts for different molecules (other than hexanal, which is the molecule used in these examples).

**Note**: To use these scripts, LAMMPS must be installed with the FEP package.

### To modify the equilibration scripts:
+ The **create_box** command needs to be modified to include the required number of atom types, bonds, and angles. Generally, you will need 3 + N atom types, 1 + B bond types, and 1 + A angle types, where N, B, and A are the numbers of atom, bond, and angle types in the fragrance molecule respectively. 
+ The **read_data** command needs to be modified. Replace hexanal_data.dat with {chemical_name}_data.dat, and hexanal_solvent_data.dat with {chemical_name}_solvent_data.dat. If you are changing the length of the polyurea chains, PU_20.dat will need to be changed to PU_{number_of_beads_per_chain}.dat
+ If you are changing the simulation box size, **fix 2** needs to be modified accordingly. The number of water molecules must also be modified in line 16 (create_atoms 	1 ...) to give the correct density of water for the new box size.
+ To generate additional trials, we varied the seed used in the **velocity** command.

### To modify the production scripts:
+ The **fix ADAPT** and **comptute fep** commands both need to be modified to include the relevant bead types. For water simulations, this will be "1 4*(3+N)", for polyurea, this will be "2*3 4*(3+N)", and for self, this will be "4*(3+N) (4+N)*(3+2N)" with N the number of bead types in the fragrance molecule.
