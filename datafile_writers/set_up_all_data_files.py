#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:12:13 2023

@author: skronen
"""

import numpy as np
import os
from martini_mapping_from_sdf import write_diffusing_molecule_data_file
from write_nonbond import write_molecule_nonbond_list, write_nonbond_list_equil
from write_pu_dat_file import write_PU_data_file
from util_functions import compute_number_of_molecs
from add_soft_potential import write_backward_difference_data_file, write_soft_nonbond_datafiles, add_additional_atypes_for_solvent


def write_data_files(chemical, grouping, mapping, len_box, save_dir, sim_type, 
                 len_box_production, num_chains, num_monomers, solvent,
                 molec_density, molec_mass, ti_pairlist, d_lambda):
    try:
        os.mkdir(save_dir)
    except:
        pass
    write_diffusing_molecule_data_file(chemical, grouping, mapping, 5, len_box, save_dir, plot_result=False)
    
    write_molecule_nonbond_list(f'Nonbonded_{chemical}.xlsx', mapping, save_dir)
    write_nonbond_list_equil(f'Nonbonded_{chemical}.xlsx', mapping, save_dir)
    
    if solvent.lower() =='water':
        molec_density = 0.997
        molec_mass = 18.02
        
        num_chains = 0
        write_PU_data_file(num_chains, num_monomers, len_box, save_dir)
        num_waters = compute_number_of_molecs(molec_density, len_box_production, molec_mass)
        print(f'{int(num_waters/4)} water beads needed for the correct density of water')
         
    elif solvent.lower() == 'pu':
        write_PU_data_file(num_chains, num_monomers, len_box, save_dir)
        num_molecs = 0
        write_diffusing_molecule_data_file(chemical, grouping, mapping, num_molecs, len_box, save_dir, plot_result=False, solvent = True)
    
    
    elif solvent.lower() == 'self':
        num_chains = 0
        write_PU_data_file(num_chains, num_monomers, len_box, save_dir)
        num_molecs = compute_number_of_molecs(molec_density, len_box_production, molec_mass) - 1
        write_diffusing_molecule_data_file(chemical, grouping, mapping, num_molecs, len_box, save_dir, plot_result=False, solvent = True)
    
    
    if sim_type == 'ti':
        #add additional atom types to the nonbond for solvent self
        add_additional_atypes_for_solvent(f'{save_dir}/nonbond_equil.txt')
        add_additional_atypes_for_solvent(f'{save_dir}/nonbond_LJ.txt')
        
        write_soft_nonbond_datafiles(f'{save_dir}/nonbond_equil_additional.txt')
        write_soft_nonbond_datafiles(f'{save_dir}/nonbond_LJ_additional.txt')
        
        write_backward_difference_data_file(f'{save_dir}/nonbond_LJ_additional_soft.txt', ti_pairlist, d_lambda)
            
"""
These two examples demonstrate how the code can be organized to write the data files 
required for the simulation.



"""

def codeine_diffusion_example():
    #Example of MARTINI model description for codeine (obtained using visualize_sdfs.py)
    chemical = 'codeine'
    grouping = [[4, 11], [2, 3],   [9,8],[5, 6, 15], [20, 17],[32, 35, 38], [16, 7, 31], [19,22], [23, 28, 27] ]
    mapping = ['TN2a', 'TC5',  'TC5',  'TN4a',  'TP1',   'SN1',     'SC3',  'TC4',  'SC3']
    
    
    '''General inputs '''
    len_box = 200 # length of initial simulation box before compression
    save_dir = 'codeine_diffusion_test' #where to save the data files
    len_box_production = 100 # length of simulation box after compression (used to calculate number of solvent molecules required for TI simulations)
    sim_type = 'diffusion' #'diffusion' or 'ti'
    
    '''Polyurea Info''' #Change the number of chains/monomers as desired
    num_chains = 442
    num_monomers = 20 #actually total number of beads (ie num_monomers = 20 -> 20 beads = 10 monomers here)
    
    '''Solvent Info'''
    solvent = 'PU'  #'water', 'PU' (polyurea), or 'self' (more molecules of the same type)
    
    #change these to the desired values for your FM:
    molec_density = 1.32 #g/cm3 codeine
    molec_mass = 299.364 #g/mol codeine
    
    '''TI Params''' #only need to be specified if sim_type == 'ti'
    #A list of the pairs that you are performing ti on. 
    #Since this is a diffusion simulation, these variables are ignored
    ti_pairlist = [[i,j] for i in range(4, 13) for j in range(2, 4)]
    d_lambda = 0.002
    
    write_data_files(chemical, grouping, mapping, len_box, save_dir, sim_type, 
                     len_box_production, num_chains, num_monomers, solvent,
                     molec_density, molec_mass, ti_pairlist, d_lambda)
    
def hexanal_ti_example():
    #Example of MARTINI model description for codeine (obtained using visualize_sdfs.py)
    chemical = 'hexanal'
    grouping = [[1,2,8,10], [14,16,20]]
    mapping = ['C1', 'SN6a']
    
    '''General inputs '''
    len_box = 200 # length of initial simulation box before compression
    save_dir = 'hexanal_ti_self' #where to save the data files
    len_box_production = 100 # length of simulation box after compression (used to calculate number of solvent molecules required for TI simulations)
    sim_type = 'ti' #'diffusion' or 'ti'
    
    '''Polyurea Info'''  #Change the number of chains/monomers as desired
    num_chains = 442
    num_monomers = 20 #actually total number of beads (ie num_monomers = 20 -> 20 beads = 10 monomers here)
    
    '''Solvent Info'''
    solvent = 'self' #'water', 'PU' (polyurea), or 'self' (more molecules of the same type)
    molec_density = 0.814 #g/cm3 hexanal
    molec_mass = 100.16 #g/mol hexanal
    
    '''TI Params'''
    #A list of the pairs that you are performing ti on. In this example, since solvent=='self', we want this list to include all 
    #pairs between the reference hexanal molecule (types 4-5) and the solvent hexanal molecules (type 6-7)
    # if you were doing ti in water, this would be [[1,j] for j in range(4, 6)], and for ti in pu, this would be [[i,j] for i in range(4, 6) for j in range(2, 4)]
    ti_pairlist = [[i,j] for i in range(4, 6) for j in range(6, 8)]
    d_lambda = 0.002 #finite difference size
    
    write_data_files(chemical, grouping, mapping, len_box, save_dir, sim_type, 
                     len_box_production, num_chains, num_monomers, solvent,
                     molec_density, molec_mass, ti_pairlist, d_lambda)
    
if __name__ == '__main__':
    codeine_diffusion_example()
    hexanal_ti_example()











