These scripts write the data files (including initial bead locations, bonds, and angles, bonded and nonbonded interaction coefficients, and bead masses) that are read by LAMMPS for the simulations completed in this project.

NOTE: we did not have time to automate the writing of the LAMMPS run scripts, so those are not generated via this code, though it would not be very difficult to do so. The steps to manually adjust the run scripts are described in their respective "example_{simulationtype}_simulation_scripts" directories.

The primary script that needs to be edited in order to write the data files is the "set_up_all_data_files.py" script (though additional changes will need to be made to the other files if future groups decide to change the simulation model, simulate different/additional molecules, etc.). 

**Generally the workflow for generating the data files for different fragrance molecules (FM) works as follows:**
1. Some manual work needs to be done first. This involves generating a structure data file (sdf) for the FM of interest. We used the Avogadro software (cited in our paper) to generate chemcial structures for codeine and hexanal, however other online resourses such as PubChem have sdf files available for download that can be used as well. The sdf file can be visualized by running the visualize_sdfs.py python script, which will allow you to confirm the sdf file appropriately describes the FM of interest.

2. The martini mapping needs to be completed. Once you have identified the Martini bead types and their associated real atoms, a list of the atom numbers (visualized using the visualize_sdfs.py script) and their associated bead types will need to be passed to the set_up_all_data_files.py script when you are ready to write your data file.

3. The nonbonded LJ parameters defined from the Martini model are read in from an excel spreadsheet. This format was chosen since it allowed all group members to collaborate in filling out the parameter table in Google Drive. These files should be named Nonbonded_{chemical_name}.xlsx and formatted as the examples we have provided for hexanal and codeine for the provided scripts to work correctly.

4. Lastly, running the set_up_all_data_files.py script as described by the comments in the script will create all the required data files needed for the simulation. 

The data files are always written such that water is bead type 1, polyurea is composed of bead types 2 and 3, and the FM of interest is bead types 4 to (4+N), with N the number of beads in the FM Martini mapping. For the 'self' solvent TI simulations, the bead types for the solvent FM molecules range from (4+N) to (4+2N)
