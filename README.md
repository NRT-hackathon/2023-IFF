# 2023-IFF

This repository is a home to the code written for the International Flavors and Fragrances (IFF) project in the 2023 University of Delaware Spring Hackathon (CHEG867) class. This project has involved running molecular simulations of different fragrance molecules in order to predict the controlled release of the fragrance molecules from polymer microcapsules (see https://sites.udel.edu/midas-nrt/graduate-trainee-timeline/courses-for-nrt-trainees/nrt-hackathon-course/ for additional information and final report). This repository contains the codes used to set up and analyze the simulations, as well as the solver for the release profile that uses the data collected from the simulations.

The repository is organized into the following groups of code, and each directory contains the additional required information needed to use the code:
+ **datafile_writers**: Python scripts used to write the data files (including initial bead locations, bonds, and angles, bonded and nonbonded interaction coefficients, and bead masses) that are read by LAMMPS for the simulations completed in this project
+ **example_diffusion_simulation_scripts**: Examples of the LAMMPS scripts used to run the diffusion simulations
+ **example_partition_simulation_scripts**: Examples of the LAMMPS scripts used to run the partition simulations
+ **release_solver**: Python scripts used to solve the macroscale diffusion model, ultimately providing the release profile as a function of time
+ **simulation_analysis**: Python scripts used to analyze simulation to obtain the diffusion coefficient from diffusion simulations and changes in free energy from partition simulations

