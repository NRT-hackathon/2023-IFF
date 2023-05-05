#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 12:18:26 2023

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import networkx as nx
from martini_write_dat_file import write_dat_file
from write_nonbond import write_molecule_nonbond_list, write_nonbond_list_equil

def get_all_34_node_sets(G):
    sets_3 = []
    sets_4 = []
    for node in G.nodes():
        for neigh1 in G.neighbors(node):
            for neigh2 in G.neighbors(neigh1):
                if node != neigh2:
                    if not ((neigh2, neigh1, node)) in sets_3:
                        sets_3.append((node, neigh1, neigh2))
                for neigh3 in G.neighbors(neigh2):
                    if node != neigh3 and node != neigh2:
                        sets_4.append((node, neigh1, neigh2, neigh3))
    return sets_3, sets_4

def visualize_sdf_to_do_martini_mapping(chemical, sdf_file = None):    
    file = f'SDF_files/{chemical}.sdf'
    if sdf_file:
        file = sdf_file
        
    df = pd.read_csv(file, header = None)

    G = nx.Graph()
    edges = []
    ct = 1
    for string in df[0]:
        l = string.split()
        p = []
        if "C" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "C")
            ct +=1
                
        elif "H" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "H")
            ct +=1
        
        elif "O" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "O")
            ct +=1
            
        elif "N" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "N")
            ct +=1
             
        elif len(l) == 7: 
            edges.append([int(l[0]), int(l[1])])

    G.add_edges_from(edges)

    ax = plt.figure().add_subplot(projection = '3d')
    heavy_atoms = []

    for i in range(1, len(G.nodes) +1):
        if G.nodes[i]['elem'] != 'H':
            heavy_atoms.append(i)
            pos = G.nodes[i]['pos']
            elem =  G.nodes[i]['elem']
            if elem == 'C':
                ax.scatter(pos[0], pos[1], pos[2], label = str(i), color = 'k')
                ax.text(pos[0], pos[1], pos[2],  '%s' % (str(i)), size=10, zorder=1, color = 'k') 
            elif elem == 'N':
                ax.scatter(pos[0], pos[1], pos[2], label = str(i), color = 'b')
                ax.text(pos[0], pos[1], pos[2],  '%s' % (str(i)), size=10, zorder=1, color = 'k') 
            elif elem == 'O':
                ax.scatter(pos[0], pos[1], pos[2], label = str(i), color = 'r')
                ax.text(pos[0], pos[1], pos[2],  '%s' % (str(i)), size=10, zorder=1, color = 'k') 
            else:
                ax.scatter(pos[0], pos[1], pos[2], label = str(i), color = 'purple')
                ax.text(pos[0], pos[1], pos[2],  '%s' % (str(i)), size=10, zorder=1, color = 'r') 
        else:
            pos = G.nodes[i]['pos']
            ax.scatter(pos[0], pos[1], pos[2], color = 'gray') 
            
    for edge in G.edges():
        p1 = G.nodes()[edge[0]]['pos']
        p2 = G.nodes()[edge[1]]['pos']
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]],[p1[2], p2[2]], color = 'darkorange')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
def write_diffusing_molecule_data_file(chemical, grouping, martini_mapping, num_molecs, len_box, save_dir, plot_result = True, solvent = False):
    file = f'SDF_files/{chemical}.sdf'
    df = pd.read_csv(file, header = None)
    
    G = nx.Graph()
    edges = []
    ct = 1
    for string in df[0]:
        l = string.split()
        p = []
        if "C" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "C")
            ct +=1
                
        elif "H" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "H")
            ct +=1
        
        elif "O" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "O")
            ct +=1
            
        elif "N" in l:
            for i in range(3):
                p.append(float(l[i]))
            G.add_node(ct, pos = p, elem = "N")
            ct +=1
             
        elif len(l) == 7: 
            edges.append([int(l[0]), int(l[1])])
    
    G.add_edges_from(edges)
    
    heavy_atoms = []
    
    for i in range(1, len(G.nodes) +1):
        if G.nodes[i]['elem'] != 'H':
            heavy_atoms.append(i)
    
    removed_edges = []
    for edge in list(G.edges()):
        group1 = None
        group2 = None
        cond1 = G.nodes[edge[0]]['elem'] == 'H'
        cond2 = G.nodes[edge[1]]['elem'] == 'H'
        if cond1 or cond2: 
            continue
        
        for i,g in enumerate(grouping):
            if edge[0] in g:
                group1 = i
            if edge[1] in g:
                group2 = i
        if group1==None: raise ValueError('check grouping labels')
        if group2==None: raise ValueError('check grouping labels')
        
        if not group1 == group2:
            G.remove_edge(edge[0], edge[1])
            removed_edges.append([group1, group2])
            
    removed_edges = np.unique(removed_edges, axis = 0)
    if len(removed_edges.shape) ==1:
        removed_edges = removed_edges[np.newaxis,:]
            
    cg_beads = list(nx.connected_components(G))
    
    #resort so they are in the same order as grouping
    cg_beads_order = []
    for atoms in cg_beads:
        ct = 0
        for atom in atoms:
            if atom in heavy_atoms and ct == 0:
                for i,gp in enumerate(grouping):
                    if atom in gp:
                        cg_beads_order.append(i)
                        ct +=1
    
    cg_beads = [x for _,x in sorted(zip(cg_beads_order, cg_beads))]
    
    
    G_martini = nx.from_edgelist(removed_edges)
    sets_3, sets_4 = get_all_34_node_sets(G_martini)
    
    mass_dict = {
        'H': 1.01,
        'C': 12.01,
        'O': 16.0,
        'N': 14.01}
    
    cg_pos = []
    cg_mass = []
    for bead in cg_beads:
        mass = 0
        pos = []
        for atom in bead:
            el = G.nodes[atom]['elem']
            m = mass_dict[el]
            mass +=m
            p = G.nodes[atom]['pos']
            pos.append(p)
        pos = np.array(pos)
        cg_pos.append(np.mean(pos, 0))
        cg_mass.append(np.round(mass,2))
    
    if plot_result:
        ax2 = plt.figure().add_subplot(projection = '3d')
        
        cg_pos = np.round(cg_pos, 4)
        for i,p in enumerate(cg_pos):
           # ax.scatter(p[0], p[1], p[2], s = 30, c = 'purple')
            ax2.scatter(p[0], p[1], p[2], s = 30, c = 'purple')
            #ax2.text(p[0], p[1], p[2],  '%s' % (str(i)), size=10, zorder=1, color = 'k') 
                    
        for edge in G_martini.edges():
            p1 = cg_pos[edge[0]]
            p2 = cg_pos[edge[1]]
            ax2.plot([p1[0], p2[0]], [p1[1], p2[1]],[p1[2], p2[2]], color = 'darkorange')
        
    bond_dists = []
    for b in removed_edges:
        p1 = cg_pos[b[0]]
        p2 = cg_pos[b[1]]
        d = np.linalg.norm(p2- p1)
        bond_dists.append(d)
    
    angles = []
    for s in sets_3:
        p1 = cg_pos[s[0]]
        p2 = cg_pos[s[1]]
        p3 = cg_pos[s[2]]
        v1 = p1 - p2
        v2 = p3 - p2
        theta = np.arccos(np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2))
        angles.append(theta/np.pi*180)
    
    nonbond_excel = f'{chemical}_Nonbonded.xlsx'
    write_dat_file(chemical, cg_pos, bond_dists, angles, cg_mass, removed_edges, sets_3, num_molecs, len_box, nonbond_excel, save_dir, solvent)




# chemical = 'codeine'
# # grouping = [[2,3,4,6], [1,5,7]] #hexanal
# # mapping = ['C1', 'SN6a']
# # save_dir = 'Hexanal_Diffusion'



# grouping = [[3,22],[20,21],[19,13],[1,16,8],[2,14],[18,4,12],[5,10,9],[17,15],[6,7,11]] #codeine
# mapping = ['TN2a', 'TC5',  'TC5',  'TN4a',  'TP1',   'SN1',     'SC3',  'TC4',  'SC3'] #codeine
# save_dir = 'test_save_dir'

# write_diffusing_molecule_data_file('codeine', mapping, 1, 100, save_dir)

# write_molecule_nonbond_list(f'Nonbonded_{chemical}.xlsx', mapping, save_dir)
# write_nonbond_list_equil(f'Nonbonded_{chemical}.xlsx', mapping, save_dir)



