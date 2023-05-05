#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:56:59 2023

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('/home/skronen/Documents/Hackathon')
from util_functions import compute_fdti

#this needs to be changed if you use a different lambda step for the finite difference
d_lambda = 0.002

ti_data = []
ti_results = []

for file in os.listdir(os.getcwd()):
    if file.endswith('.fep'):
        data = np.loadtxt(file)   
        plt.plot(data[:,0], data[:,1])
        ti_data.append(data[:,1])
        result = compute_fdti(file, d_lambda)
        ti_results.append(result)
        
timesteps = data[:,0]
ti_data = np.array(ti_data)

lambs = np.linspace(0,1,len(timesteps))

derivs = ti_data/d_lambda
means = np.mean(derivs, 0)*4.18
stds = np.std(derivs, 0)*4.18

plt.figure()  
plt.errorbar(lambs, means, yerr  = stds, color = 'k')

colors = ['red', 'blue']
for i in range(len(lambs)-1):
    c_ind = i%2
    plt.fill_between(lambs[i:i+2], 0, means[i:i+2], color = colors[c_ind], alpha = 0.1)

plt.figure()  
plt.errorbar(lambs, means, yerr  = stds, color = 'r')

plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\dfrac{dU}{d\lambda}$ (kJ/mol)')

ti_results = np.array(ti_results)
ti_kJ = ti_results * 4.18
ti_mean = np.mean(ti_kJ)
ti_std = np.std(ti_kJ)

print(f'delG = {ti_mean:0.2f} ({ti_std:0.2f})')