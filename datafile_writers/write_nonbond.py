#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:06:56 2023

@author: skronen
"""

import numpy as np
import pandas as pd

def write_nonbond_coeffs(excel_file):
    sigmas = pd.read_excel(excel_file, sheet_name='sigmas', header = 1).to_numpy()
    epsilons = pd.read_excel(excel_file, sheet_name='epsilons', header = 1).to_numpy()
    
    num_beadtypes = len(sigmas)
    
    nonbond_list = []

    for i in range(num_beadtypes):
        for j in range(i, num_beadtypes):
            eps = np.round(epsilons[i][j] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas[i][j] * 10,3) #convert from nm to angstroms
            nonbond_list.append(f'{i+1} {j+1} {eps} {sig}')
            nonbond_list.append('\n')
            
    return nonbond_list

def write_nonbond_w_pu(excel_file, WCA = False):
    sigmas = pd.read_excel(excel_file, sheet_name='sigmas', header = 1).to_numpy()
    epsilons = pd.read_excel(excel_file, sheet_name='epsilons', header = 1).to_numpy()
    
    nonbond_list = []

    for i in range(3):
        for j in range(i, 3):
            eps = np.round(epsilons[i][j] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas[i][j] * 10,3) #convert from nm to angstroms
            if WCA:
                nonbond_list.append(f'{i+1} {j+1} {eps} {sig} {np.round(sig*2**(1/6),3)}')
            else: 
                nonbond_list.append(f'{i+1} {j+1} {eps} {sig}')
            nonbond_list.append('\n')
            
    return nonbond_list

def write_molecule_nonbond_list(excel_file, mapping_list, save_dir):

    sigmas = pd.read_excel(excel_file, sheet_name='sigmas', header = 1)
    sigmas_np = sigmas.to_numpy()
    epsilons = pd.read_excel(excel_file, sheet_name='epsilons', header = 1)
    epsilons_np = epsilons.to_numpy()
    
    beads = sigmas.columns

    num_beadtypes = len(mapping_list)
    
    w_PU_nonbond_list = write_nonbond_w_pu(excel_file)
    
    nonbond_list = []
    
    #add water/PU interactions
    for j,bead_type in enumerate(mapping_list):
        for i in range(3):
            bead1_ind = np.where(beads == bead_type)[0][0]
            bead2_ind = i
            if not bead1_ind > bead2_ind:
                bead1_ind, bead2_ind = bead2_ind, bead1_ind

            eps = np.round(epsilons_np[bead2_ind][bead1_ind] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas_np[bead2_ind][bead1_ind] * 10,3) #convert from nm to angstroms
            nonbond_list.append(f'{i+1} {j+4} {eps} {sig}')
            nonbond_list.append('\n')

        
    #add self interactions
    for i,bead1 in enumerate(mapping_list):
        for j,bead2 in enumerate(mapping_list[i:]):
            j += i
            bead1_ind = np.where(beads == bead1)[0][0]
            bead2_ind = np.where(beads == bead2)[0][0]
            if not bead1_ind < bead2_ind:
                bead1_ind, bead2_ind = bead2_ind, bead1_ind
            eps = np.round(epsilons_np[bead1_ind][bead2_ind] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas_np[bead1_ind][bead2_ind] * 10,3) #convert from nm to angstroms
            nonbond_list.append(f'{i+4} {j+4} {eps} {sig}')
            nonbond_list.append('\n')
    
    with open(f'{save_dir}/nonbond_LJ.txt', 'w') as f:
        f.write('\n')
        f.write('0 atoms \n')
        f.write(f'{num_beadtypes + 3} atom types \n')
        f.write('\n')
            
        f.write('PairIJ Coeffs\n\n')
        for string in w_PU_nonbond_list:
            f.write(string)
        
        for string in nonbond_list:
            f.write(string)
        f.write('\n')
        
def write_nonbond_list_equil(excel_file, mapping_list, save_dir):

    sigmas = pd.read_excel(excel_file, sheet_name='sigmas', header = 1)
    sigmas_np = sigmas.to_numpy()
    epsilons = pd.read_excel(excel_file, sheet_name='epsilons', header = 1)
    epsilons_np = epsilons.to_numpy()
    
    beads = sigmas.columns

    num_beadtypes = len(mapping_list)
    
    w_PU_nonbond_list = write_nonbond_w_pu(excel_file, WCA = True)
    
    nonbond_list = []
    #add water/PU interactions
    for j,bead_type in enumerate(mapping_list):
        for i in range(3):
            bead1_ind = np.where(beads == bead_type)[0][0]
            bead2_ind = i
            if not bead1_ind > bead2_ind:
                bead1_ind, bead2_ind = bead2_ind, bead1_ind

            eps = np.round(epsilons_np[bead2_ind][bead1_ind] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas_np[bead2_ind][bead1_ind] * 10,3) #convert from nm to angstroms
            nonbond_list.append(f'{i+1} {j+4} {eps} {sig} {np.round(sig*2**(1/6),3)}')
            nonbond_list.append('\n')

        
    #add self interactions
    for i,bead1 in enumerate(mapping_list):
        for j,bead2 in enumerate(mapping_list[i:]):
            j += i
            bead1_ind = np.where(beads == bead1)[0][0]
            bead2_ind = np.where(beads == bead2)[0][0]
            if not bead1_ind < bead2_ind:
                bead1_ind, bead2_ind = bead2_ind, bead1_ind
                
            eps = np.round(epsilons_np[bead1_ind][bead2_ind] /4.18, 3) #convert from kJ to kcals
            sig = np.round(sigmas_np[bead1_ind][bead2_ind] * 10,3) #convert from nm to angstroms
            nonbond_list.append(f'{i+4} {j+4} {eps} {sig} {np.round(sig*2**(1/6),3)}')
            nonbond_list.append('\n')
            

    
    with open(f'{save_dir}/nonbond_equil.txt', 'w') as f:
        f.write('\n')
        f.write('0 atoms \n')
        f.write(f'{num_beadtypes + 3} atom types \n')
        f.write('\n')
            
        f.write('PairIJ Coeffs\n\n')
        for string in w_PU_nonbond_list:
            f.write(string)
        
        for string in nonbond_list:
            f.write(string)
        f.write('\n')