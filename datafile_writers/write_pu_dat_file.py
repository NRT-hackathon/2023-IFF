#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:56:31 2023

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from write_nonbond import write_nonbond_w_pu

def snake2d(N, L, spacing):
    L = int((L-1)/spacing)
    num_cols = N//L + 1* (N%L>0)
    num_rows = L-1
    Nf = num_cols * num_rows
    arr = np.reshape(np.linspace(0, Nf-1, Nf), (num_cols, num_rows))
    out = arr.copy()
    out[1::2] = arr[1::2,::-1]
    out = out.ravel()
    x = np.linspace(1, num_cols, num_cols)
    y = np.linspace(1, num_rows, num_rows)
    pos = list(product(x,y))
    pos = [x for _, x in sorted(zip(out,pos))]
    pos = np.array(pos)
    pos = pos[:N]
    z = np.expand_dims(np.zeros(N), 1)
    pos = np.append(pos, z, axis = 1) 
    pos = pos* spacing
    return pos   

def write_PU_data_file(num_chains, num_monomers, len_box, save_dir):
    for L_num, L_val in enumerate([20]):
        filename = f"{save_dir}/PU_" + str(L_val) + ".dat"
        N = num_monomers #number of beads
        M = num_chains #number of chains
        bond_length = 4.3
        
        dat_atoms = np.zeros([M* N, 6], dtype = object)  #id molecid atomtype  x y z
        dat_bonds= np.zeros([M*(N-1), 4], dtype= object) #id type atom1 atom2
        dat_angles = np.zeros([M*(N-2), 5], dtype = object) #id type atom1 atom2 atom3
        one_poly = snake2d(N, len_box, bond_length) #initialize reference polymer that will be rotated and translated
          
        space = one_poly[-1,0] - one_poly[0,0] + bond_length
        grid_size = int(len_box /space)
        x_spacing = len_box/grid_size
        y_spacing = bond_length
        
        for i in range(M):         
            new_poly = one_poly.copy()
            new_poly[:,0] += x_spacing * (i%grid_size)
            new_poly[:,2] += y_spacing * (i//grid_size)
            dat_atoms[N*i: N*(i+1), 3:6] = new_poly
          
        """
        ATOMS
        """
        for i in range(M*N):
            dat_atoms[i,0] = (i+1)
            dat_atoms[i,1] = int(i/N)+1
            if i%2 ==0:
                dat_atoms[i,2] = 1
            else:
                dat_atoms[i,2] = 2
        
        """
        BONDS
        """
        for i in range(M*(N-1)):
            dat_bonds[i,0] = i+1
            dat_bonds[i,1] = 1
        
        count=0
        for i in range(M):
            for j in range(N-1):
                dat_bonds[count,2] = (j+1 + N*i)
                dat_bonds[count,3] = (j+2 + N*i)
                count +=1
                
        """
        ANGLES
        """
        if N>2:
            for i in range(M*(N-2)):
                dat_angles[i,0] = i+1
                dat_angles[i,1] = 1
            
            count=0
            for i in range(M):
                for j in range(N-2):
                    dat_angles[count,2] = (j+1 + N*i)
                    dat_angles[count,3] = (j+2 + N*i)
                    dat_angles[count,4] = (j+3 + N*i)
                    count +=1
    
    
        if N<=2:  
            raise ValueError('Pick number of monomers greater than 2')
            # """
            # WRITE TO TEXT FILE
            # """
            # top = ['#Data File for Basic Polymer Chains',
            #        '{} atoms'.format(M*N),
            #        '{} bonds'.format((M*(N-1))),
            #        '',
            #        '{} atom types'.format(2),
            #        '{} bond types'.format(1),
            #        '',
            #        '-{} {} xlo xhi'.format(len_box/2,len_box/2),
            #        '-{} {} ylo yhi'.format(len_box/2,len_box/2),
            #        '-{} {} zlo zhi'.format(len_box/2,len_box/2),
            #        '',
            #        'Masses',
            #        '#id mass',
            #        '1 56.0',
            #        '2 58.0',
            #        '']       
            
            # #Write top to file
            # textfile = open(filename, "w")
            # for element in top:
            #     textfile.write(element + "\n")
            
            # textfile.write('Bond Coeffs\n\n')
            # textfile.write('1 299 4.3')
            # textfile.write('\n\n')
    
            # textfile.write("Atoms" + "\n" + "#id molectype atomtype x y z" + "\n")
            # for element in dat_atoms:
            #     textfile.write(str(element)[1:-1] + "\n")
            
            # textfile.write("\n" + "Bonds" + "\n" + "#id type atom1 atom2" + "\n")
            # for element in dat_bonds:
            #     textfile.write(str(element)[1:-1] + "\n")
            
            # textfile.close()
            
        
        else:
            """
            WRITE TO TEXT FILE
            """
            
                
            top = ['#Data File for Polyurea Chains',
                   '{} atoms'.format(M*N),
                   '{} bonds'.format((M*(N-1))),
                   '{} angles'.format((M*(N-2))),
                   '',
                   '{} atom types'.format(2),
                   '{} bond types'.format(1),
                   '{} angle types'.format(1),
                   '',
                   '-{} {} xlo xhi'.format(len_box/2,len_box/2),
                   '-{} {} ylo yhi'.format(len_box/2,len_box/2),
                   '-{} {} zlo zhi'.format(len_box/2,len_box/2),
                   '',
                   'Masses',
                   '#id mass',
                   '1 56.0',
                   '2 58.0',
                   '']       
            
            
            
            #Write top to file
            textfile = open(filename, "w")
            for element in top:
                textfile.write(element + "\n")
            
            textfile.write('Bond Coeffs\n\n')
            textfile.write('1 299 4.3')
            textfile.write('\n\n')
    
            textfile.write('Angle Coeffs\n\n')
            textfile.write('1 224.9 142.12')
            #textfile.write('1 36 180')
            textfile.write('\n\n')
            
            if num_chains >0:
                textfile.write("Atoms" + "\n" + "#id molectype atomtype x y z" + "\n")
                for element in dat_atoms:
                    textfile.write(str(element)[1:-1] + "\n")
                
                textfile.write("\n" + "Bonds" + "\n" + "#id type atom1 atom2" + "\n")
                for element in dat_bonds:
                    textfile.write(str(element)[1:-1] + "\n")
                
                textfile.write("\n" + "Angles" + "\n" + "#id type atom1 atom2 atom3" + "\n")
                for element in dat_angles:
                    textfile.write(str(element)[1:-1] + "\n")
            textfile.close()