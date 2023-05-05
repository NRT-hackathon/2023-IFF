#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:00:05 2023

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
from Dump_Reader import read_dump
import os

f = os.listdir(os.getcwd())
files = []
for file in f:
    if file.endswith('.lammpstrj'):
        files.append(file)
files.sort()

ax1 = plt.figure().add_subplot()
all_trials = []

#masses of beads in codeine, this needs to be changed for a different molecule
mass = {4:31.04, 
        5: 25.03, 
        6: 25.03, 
        7: 41.03, 
        8: 30.03, 
        9: 43.08, 
        10: 38.05, 
        11: 26.04, 
        12: 40.07}


for file in files:
    all_molecs = []
    dump = read_dump(file)
    
    len_box = dump[0][0]
    timesteps = np.array(dump[1])
    timesteps = (timesteps - np.min(timesteps)) * 3 #timesteps in fs
    xpos = dump[2]
    ypos = dump[3]
    zpos = dump[4]
    atypes = dump[5]
    molecs = dump[6]
    
    for m in np.unique(molecs):
        mask = molecs == m
        
        relevant_types = atypes[mask]
        masses = np.array([mass[typ] for typ in relevant_types])
        masses = masses/np.sum(masses)
        relevant_x = xpos.T[mask][:,np.newaxis]
        relevant_y = ypos.T[mask][:,np.newaxis]
        relevant_z = zpos.T[mask][:,np.newaxis]
        
        dt = timesteps[1] - timesteps[0]
        all_slopes = []
        for lag in np.arange(0, len(timesteps)):
            molec_pos = np.concatenate([relevant_x, relevant_y, relevant_z], axis = 1)
            molec_pos = molec_pos[:,:,lag:]
            mass_weighted_pos = masses[:,np.newaxis, np.newaxis] * molec_pos
            com = np.sum(mass_weighted_pos, axis = 0)
            normal_pos = com - com[:,0][:,np.newaxis] #normalize so initial position is the origin
            msd = np.sum(normal_pos**2, axis = 0)
            msd = msd/100 #convert from A2 to nm2
            
            slopes = []
            for m in msd[1:]:
                slopes.append((m - msd[0]))
            all_slopes.append(slopes)    
                
            if lag ==0:
                msd_full = msd
            
        msd = msd_full
        
        all_slopes_arr = []
        for i,s in enumerate(all_slopes[:-1]):
            for j in range(i):
                s.append(np.nan)
            all_slopes_arr.append(s)
        
        all_slopes_arr = np.array(all_slopes_arr)
        
        slopes = np.nanmean(all_slopes_arr, 0)
        errs = np.nanstd(all_slopes_arr, 0)
        
        all_molecs.append(slopes)
        
        ax1.loglog(timesteps[:-1], slopes)
        max_ts = np.log10(np.max(timesteps))
        x = np.logspace(6, max_ts, 100)
        y = x / 2e5
        ax1.plot(x, y, ls = '--')
        plt.ylabel(r'msd ($nm^2$)')
        plt.xlabel('time (fs)')
    all_trials.append(all_molecs)

all_trials = np.array(all_trials)
mean = np.mean(all_trials, 1)
std = np.std(all_trials, 1)
times =timesteps[1:]
log_ts = np.log10(times)

i = 0
D_vals = []
for m,s in zip(mean, std):
    plt.figure()
    plt.loglog(times, m, color = 'r')
    plt.fill_between(times, m - s, m + s,color = 'r', alpha = 0.5)
    x = np.logspace(np.log10(times[0]), max_ts, 100)
    y = x / 3.5e5 #this is just a random scaling factor
    plt.plot(x, y, ls = '--', color = 'k')
    plt.ylabel(r'msd ($nm^2$)')
    plt.xlabel('time (fs)')
        
    log_slopes = np.log10(all_trials[i][0])
    
    linear_mask = np.logical_and(times>0, times<1.5e7)
    m,b = np.polyfit(times[linear_mask], m[linear_mask], 1)
    
    m_si = m/1000
    i+=1
    D_vals.append(m_si)
    
D_mean = np.mean(D_vals)
D_std = np.std(D_vals)
exponent = int(np.floor(np.log10(D_mean)))
D_mean /= 10**exponent
D_std /= 10**exponent


print(f'D is {D_mean:0.3f} ({D_std:0.3f}) e{exponent} m^2/s')    
