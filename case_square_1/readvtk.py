#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 12:14:41 2025

@author: lsicard
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import pyvista as pv
import os
import re

def natural_key(string):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string)]

test_location = 11036

folder_path = "/home/lsicard/peka5/Peka2D_v5-1/case_square_1"

os.chdir("/home/lsicard/peka5/Peka2D_v5-1/case_square_1")

vtk_files = [f for f in os.listdir(folder_path) if f.endswith('.vtk')]
vtk_files.sort(key=natural_key)

nb_solute = 6


##Get first vqlue of phi
list_sol = []
import pyvista as pv

mesh = pv.read(vtk_files[0])
Nb_point = mesh.n_points
Nb_cells = mesh.n_cells
points = mesh.points
cells = mesh.cells
data = mesh.cell_data.keys()
data_solute = data[len(data)-nb_solute:]

for i in data_solute:
    
    if i in mesh.cell_data:
        phi = mesh.cell_data[i]
    
    list_sol.append(phi)


### get all the values of phi
for file in vtk_files[1:]:

    mesh = pv.read(file)
    
    Nb_point = mesh.n_points
    Nb_cells = mesh.n_cells
    
    points = mesh.points
    
    cells = mesh.cells
    
    data = mesh.cell_data.keys()
    data_solute = data[len(data)-nb_solute:]
    
    counter = 0
    for i in data_solute:
        
        if i in mesh.cell_data:
            phi = mesh.cell_data[i]
        
        list_sol[counter] = np.append(list_sol[counter],phi)
        counter+=1 
        
        
##plot one solution

##Get the values at one value

DataPerPoint = []
x = np.arange(1,7,1)
for time in range(1,len(vtk_files)+1):
    DataPerPoint = []
    for l in list_sol:
        DataPerPoint.append(l[test_location*time])
    plt.plot(DataPerPoint, x)
    plt.gca().invert_yaxis()
        

    
