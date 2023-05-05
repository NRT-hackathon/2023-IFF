#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 18:46:31 2022

@author: skronen
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

#read dump file
def read_dump(file):
    if len(sys.argv) > 1:
        file = sys.argv[1]
        
    with open(file, 'r') as f:
            text = f.read()
    lines = text.splitlines()
    
    atoms = int(lines[3].split()[0])
    ts = 10000 #if you every have a trajectory with more than this many timesteps, this will need to be changed
    xpos = np.zeros([ts, atoms])
    ypos = np.zeros([ts, atoms])
    zpos = np.zeros([ts, atoms])
    
    timesteps = []
    len_box = []
    atype= np.zeros([atoms], dtype = 'int')
    molecnum= np.zeros([atoms], dtype = 'int')

    natoms = 0
    
    count = 0
    for l_num, line in enumerate(lines):
        line = line.split()
        if (len(line) ==6) and (line[0] != "ITEM:"):  
            x = float(line[3])
            y = float(line[4])
            z = float(line[5])            
            if count<atoms:
                atype[count]= int(line[2])
                molecnum[count] = int(line[1])
            
            xpos[int(count/atoms), (count - atoms*int(count/atoms))]=x
            ypos[int(count/atoms), (count - atoms*int(count/atoms))]=y
            zpos[int(count/atoms), (count - atoms*int(count/atoms))]=z
            count += 1
            
        elif (len(line) ==2) and (line[1] == "TIMESTEP"):
            timesteps.append(lines[l_num+1])
            if (len(timesteps) == 2):
                natoms = count
        elif (len(line) ==6) and (line[1] =="BOX"):
            a = lines[l_num+1].split()[1]
            len_box.append(float(a))
        else:
            continue
    xcnt = 0
    ycnt = 0
    zcnt = 0
    if np.all(xpos[0,:]==0):
        xcnt =1
    if np.all(ypos[0,:]==0):
        ycnt = 1
    if np.all(zpos[0,:]==0):
        zcnt = 1
        
    xpos = xpos[~np.all(xpos == 0, axis=1)]
    ypos = ypos[~np.all(ypos == 0, axis=1)]
    zpos = zpos[~np.all(zpos == 0, axis=1)]
    timesteps = [int(ts) for ts in timesteps]
    if xcnt ==1:
        xpos = np.vstack((np.zeros(atoms), xpos))
    if ycnt ==1:
        ypos = np.vstack((np.zeros(atoms), ypos))
    if zcnt ==1:
        zpos = np.vstack((np.zeros(atoms), zpos))
    return(len_box, timesteps, xpos, ypos, zpos, atype, molecnum, atoms)


