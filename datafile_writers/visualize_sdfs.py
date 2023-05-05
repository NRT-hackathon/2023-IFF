#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 10:23:29 2023

@author: skronen
"""

from martini_mapping_from_sdf import visualize_sdf_to_do_martini_mapping

files = ['SDF_files/hexanal.sdf', 'SDF_files/hexanal_old.sdf']

for file in files:
    visualize_sdf_to_do_martini_mapping('codeine', sdf_file = file)