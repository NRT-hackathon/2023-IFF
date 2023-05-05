#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:55:31 2023

@author: skronen
"""

import numpy as np
import sys
import os

def write_dat_file(chemical, cg_pos, bond_dists, angles, cg_mass, removed_edges, sets_3, N, len_box, nonbond_excel, save_dir, solvent = False):
    spacing = np.max(cg_pos, 0) - np.min(cg_pos, 0) + 6
    
    num_copies = (len_box /spacing).astype(int)
    layer_size = num_copies[0] * num_copies[1]
        
    num_beads = len(cg_pos)
    filename = f'{save_dir}/{chemical}_data.dat'
    
    if solvent == True:
        filename = f'{save_dir}/{chemical}_solvent_data.dat'

    with open(filename, 'w') as f:
        f.write(f'{chemical} data file\n')
        f.write(f'{N*len(cg_pos)} atoms \n')
        f.write(f'{len(cg_pos)} atom types \n')
        
        if num_beads >1:
            f.write(f'{N*len(bond_dists)} bonds \n')
            f.write(f'{len(bond_dists)} bond types \n')

        if num_beads > 2:
            f.write(f'{N*len(angles)} angles \n')
            f.write(f'{len(angles)} angle types \n')
        
        f.write('\n')
        
        f.write(f'{-len_box/2} {len_box/2} xlo xhi \n{-len_box/2} {len_box/2} ylo yhi \n{-len_box/2} {len_box/2} zlo zhi\n\n')
        
        f.write('Masses\n\n')
        for i,m in enumerate(cg_mass):
            f.write(f'{int(i+1)} {m}\n')
        f.write('\n')
        
        # f.write('PairIJ Coeffs\n\n')
        # for string in nonbond_list:
        #     f.write(string)
        # f.write('\n')
        
        if num_beads >1:
            f.write('Bond Coeffs\n\n')
            for i,b in enumerate(bond_dists):
                f.write(f'{i+1} {1250/4.18} {np.round(b,2)} \n')
            f.write('\n')

            
        if num_beads >2:
            f.write('Angle Coeffs\n\n')
            for i,a in enumerate(angles):
                f.write(f'{i+1} {50/4.18} {np.round(a,2)} \n')
            f.write('\n')
       
        #dihedrals are K =10kJ/mol 
        #not dealing with dihedrals for the time being
        if N>0:
            f.write('Atoms\n\n')
            for j in range(N):
                shift = np.array([j%num_copies[0], (j -j//layer_size*layer_size)  // num_copies[0], j//layer_size]) * spacing
                for i,p in enumerate(cg_pos):
                    f.write(f'{int(num_beads*j + i+1)} {j+1} {int(i+1)} {np.round(p[0] + shift[0],2)} {np.round(p[1] + shift[1],2)} {np.round(p[2] + shift[2],2)} \n')
            f.write('\n')  
            
            if num_beads >1:
                f.write('Bonds\n\n')
                for j in range(N): 
                    for i,b in enumerate(removed_edges):
                        f.write(f'{num_beads*j + int(i+1)} {int(i+1)} {num_beads*j + b[0]+1} {num_beads*j + b[1]+1} \n')
                f.write('\n')
    
            if num_beads >2:
                f.write('Angles\n\n')
                for j in range(N):
                    for i,b in enumerate(sets_3):
                        f.write(f'{num_beads*j + int(i+1)} {int(i+1)} {num_beads*j + b[0]+1} {num_beads*j + b[1]+1} {num_beads*j + b[2] + 1} \n')
                f.write('\n')
            
      
    