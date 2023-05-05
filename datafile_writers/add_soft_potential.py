#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:17:59 2023

@author: skronen
"""

import numpy as np

def add_additional_atypes_for_solvent(nonsoft_file):
    file = nonsoft_file
    
    original_lines = []
    additional_lines = []
    with open(file, 'r') as f:
        lines = f.read().splitlines()
        
        atypes = []
        for i,line in enumerate(lines):
            l = line.split(' ')
            if l[0].isnumeric() and l[1].isnumeric():
                atypes.append(float(l[0]))
                atypes.append(float(l[1]))
        max_atype = max(atypes)     
        
        for i,line in enumerate(lines):
            if 'atom types' in line:
                l = line.split(' ')
                orig_num = int(l[0])
                new_num = orig_num + max_atype - 3
                line = f'{int(new_num)} atom types'
                original_lines.append(line)
                continue
            original_lines.append(line)
            write = False
            l = line.split(' ')
            if l[0].isnumeric() and l[1].isnumeric():

                type1 = float(l[0])
                type2 = float(l[1])
                if type1>3 and type2>3:
                    new_type1 = type1 - 3 + max_atype
                    new_type2 = type2 - 3 + max_atype
                    l_copy = l.copy()
                    l_copy[0] = str(int(new_type1))
                    l_copy[0], l_copy[1] = l_copy[1], l_copy[0]
                    additional_lines.append(' '.join(l_copy))
                    
                    if type1 != type2:
                        l_copy = l.copy()
                        l_copy[1] = str(int(new_type2))
                        additional_lines.append(' '.join(l_copy))

                    l_copy = l.copy()
                    l_copy[0] = str(int(new_type1))
                    l_copy[1] = str(int(new_type2))
                    additional_lines.append(' '.join(l_copy))

                else:
                    if type1>3:
                        write = True
                        type1 = type1 - 3 + max_atype
                        l[0] = str(int(type1))
                    if type2>3:
                        write = True
                        type2 = type2 - 3 + max_atype
                        l[1] = str(int(type2))
                    
                    if write:
                        additional_lines.append(' '.join(l))
            
    file = file[:-4] + '_additional.txt'
    
    while original_lines[-1] == '':
        original_lines = original_lines[:-1]
        
    with open(file, 'w') as f:
        for line in original_lines:
            f.write(line)
            f.write('\n')
        for line in additional_lines:
            f.write(line)
            f.write('\n')
    
def write_soft_nonbond_datafiles(nonsoft_file):
    file = nonsoft_file
    if 'equil' in file:
        #if writing the file for equilibration, we need to add the lambda parameter
        #between the sigma and rcut
        lamlines = []
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            
            for i,line in enumerate(lines):
                l = line.split(' ')
                write_line = str()
                if len(l)>=4 and i > 5:
                    for j,item in enumerate(l):
                        write_line += item
                        write_line += ' '
                        if j == 3:
                            write_line += '1.0 '
                    lamlines.append(write_line)
                else: lamlines.append(line)
                    
                
                
        file = file[:-4] + '_soft.txt'
        
        with open(file, 'w') as f:
            for line in lamlines:
                f.write(line)
                f.write('\n')
                
    else:
        #if writing the file for production, we just need to add lambda at the end
        #of each line
        lamlines = []
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            
            for i,line in enumerate(lines):
                if len(line.split(' '))>=4 and i > 5:
                    line += ' 1.0'
                lamlines.append(line)
                
        file = file[:-4] + '_soft.txt'
        
        with open(file, 'w') as f:
            for line in lamlines:
                f.write(line)
                f.write('\n')

def write_backward_difference_data_file(data_file, ti_pairlist, d_lambda):
    with open(data_file, 'r') as f:
        lines = f.read().splitlines()
    
    write_lines = []
    for i,line in enumerate(lines):
        line = line.split(' ')
        if len(line)==5 and i > 5:
            type1 = int(line[0])
            type2 = int(line[1])
            if [type1, type2] in ti_pairlist or [type2, type1] in ti_pairlist:
                line[4] = str(d_lambda)
        write_lines.append(line)
    write_lines = [' '.join(w) for w in write_lines]

    with open(data_file[:-4] + '_final_perturb.txt', 'w') as f:
        for line in write_lines:
            f.write(line)
            f.write('\n')


