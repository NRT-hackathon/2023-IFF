# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 10:58:24 2023

@author: steph
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import time

def do_the_solving(params, t_max, n_times):
    #unpack parameters
    [r_mean, r_std, Cd_init, Vd, D, 
     Kmd, Kmr, h, Vr, rmin, rmax, rstep] = params
    
    """Generate Distribution of Capsule Sizes
    and their respective proportions accoring to 
    a normal distribution"""
    sizes = np.arange(rmin, rmax + rstep, rstep)
    #the following line is for debugging- uncomment to only look at the mean capsule size
    #sizes = np.array([r_mean]) 
    Npi = []
    for s in sizes:
        cdf_pts = [s-rstep/2, s+rstep/2]
        cdf_vals = stats.norm.cdf(cdf_pts, loc = r_mean, scale = r_std)
        Npi.append(cdf_vals[1] - cdf_vals[0])
    
    Npi = np.array(Npi)
    Ai = 4*np.pi*sizes**2
    Vi = 4/3*np.pi*sizes**3
    
    Vd_calc = np.sum(Npi*Vi)
    Npi = Npi*Vd/Vd_calc
    
    """Solve differential equations: Actually analytically derive the solution for 
    Cd(t) and just analytically solve for each t"""
    t_step = t_max/n_times #1000 timestpes
    
    times = np.arange(0, t_max, t_step)
    
    Cds = []
    
    for t in times:
        prefactor = -D*Ai/h/Vi * Kmd * Cd_init
        exponential = -D*Ai*Kmd/Vr/h * (Kmr/Kmd + Vr/Vi)
        dCdt = prefactor * np.exp(exponential * t)
        Cdi = 1/exponential * dCdt + Cd_init - prefactor/exponential #analytic solution based on Eq.2
        Cds.append(Cdi)
    
    Cds = np.array(Cds) 
    #Cds is an array of shape n_timesteps x n_capsule_radii(sizes)
    #Each row is a specific time, and each column is a specific capsule radius
    
    """Compute total release from all the capsules"""
    M_init = np.sum(Cd_init*Vi*Npi)
    delM = Cds*Vi*Npi
    delM_tot = np.sum(delM, axis = 1)
    release = 100*(M_init - delM_tot)/M_init
    return times, release

def get_bootstrap_params(mean, std):
    randvals = np.random.randn(len(mean))*std
    params = mean.copy()
    params[5:7] = np.log10(params[5:7]) 
    params += randvals
    params[5:7] = 10**params[5:7]
    return params

def bootstrap_solve(mean, std, N, t_max, n_times):
    bootstrap_results = []
    for i in range(N):
        params = get_bootstrap_params(mean, std)
        times, release = do_the_solving(params, t_max, n_times)
        bootstrap_results.append(release)
    bootstrap_results = np.array(bootstrap_results)
    times = np.array(times)
    return times, bootstrap_results
    
"""Set up parameters"""
r_mean = 129e-9 #m
r_std = 30e-9 #m

Cd_init = 358.44e3 #g/m3
Vd = 0.485e-6 #m3
D = 1.027e-19 #m2/s
Kmd = 0.11 
Kmr = 2.67
h = 2.86e-9 #m
Vr = 400e-6 #m3

rmin = 29e-9 #m    
rmax = 329e-9 # m     
rstep = 1e-8 #m

#parameters used in Scenario 1 of Muro Sune paper
params_1 = [r_mean, r_std, Cd_init, Vd, D,
                Kmd, Kmr, h, Vr, rmin, rmax, rstep]
params_2 = params_1.copy()
params_3 = params_1.copy()

#replace D/Ks with those computed from simulation (params2 = codeine, params3 = hexanal)
params_2[2] = 1200 #C_init
params_2[4] = 1.954e-9 #D
params_2[5] = 10**4.25#Kmd
params_2[6] = 10**12.9#Kmr
params_3[2] = 814 #C_init
params_3[4] = 7.16e-9 #D 
params_3[5] = 10**2.40#Kmd
params_3[6] = 10**5.70#Kmr

all_params = [params_2, params_3] #all params is the list of parameters we will loop through to solve
#These are lists of the errors we get from simulation
#std of D are absolute
#std of Kmd and Kmr are actually of log10(Kmd) and log10(Kmr)
stds_2 = [0,0,0,0,8.3e-11, 0.9, 0.73, 0 , 0 , 0, 0 , 0]
stds_3 = [0,0,0,0,1.86e-9, 0.39, 0.38, 0 , 0 , 0, 0 , 0]
all_stds = [stds_2, stds_3]

if len(all_params) != len(all_stds):
    raise ValueError('Different numbers of parameters and errors input')
    
""" Set up Plots"""
ax = plt.figure().add_subplot()
ax.set_xlabel('Time (seconds)')
ax.set_ylabel('% Release')
ax.set_title('Release Profile')

colors = ['r', 'b']
labels = ['codeine', 'hexanal']
i = 0
    
for params, std in zip(all_params, all_stds):
    tic = time.perf_counter()
    times, release = bootstrap_solve(np.array(params), np.array(std), N= 1_000, t_max = 1e-10, n_times = 100)
    toc = time.perf_counter()
    print(f'bootstrapping took {toc - tic:0.4f} seconds')
    mean_release = np.mean(release, 0)
    std_release = np.std(release, 0)
    
    #plot release
    ax.plot(times, mean_release, color =colors[i], label = labels[i]) 
    ax.fill_between(times, mean_release - std_release, mean_release + std_release, color = colors[i], alpha = 0.3)
    i+=1
ax.legend(loc = 'upper left')
ax.set_ylim([0,100])
