#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 16:41:30 2023

@author: skronen
"""

import numpy as np

def compute_number_of_molecs(molec_density, len_box, molec_mass):
    """
    Parameters
    ----------
    molec_density : float
        bulk density of molecule in g/cm3
    len_box : float
        box side length of simulation box
    molec_mass : float
        molar mass of the molecule in g/mol

    Returns
    -------
    molec_num
        The number of molecules needed to populate the simulation box

    """
    v_box = len_box**3
    unit_convert = 0.6022 #N_av = 6.022e23 / (1e24 A3/cm3)
    molec_num = v_box * molec_density * unit_convert /molec_mass
    return int(molec_num)

# N = compute_number_of_molecs(1.32, 48.354, 299.4)

# print(N)
# N = compute_number_of_molecs(1, 48.354, 72)
# print(N)

def compute_fdti(fdti_file, d_lambda):
    ti_data = np.loadtxt(fdti_file)
    
    Udiffs = ti_data[:,1]
    
    ti_sum = 0
    for i,dU in enumerate(Udiffs[:-1]):
        low = Udiffs[i]
        high = Udiffs[i+1]
        ti_sum += (high + low) / (2 * d_lambda)


    return (ti_sum/(i - 1))    # divide by i - 1 == multiply by delta

# files = [f'codeine_water_fdtiresults/res{i}.fep' for i in range(1,4)]
# for f in files:
#     res = compute_fdti(f, d_lambda=0.002)
#     print(res)

