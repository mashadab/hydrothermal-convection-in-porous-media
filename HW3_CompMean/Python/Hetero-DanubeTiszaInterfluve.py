# -*- coding: utf-8 -*-
"""
Danube Tiscza Interfluve (heterogeneous case)
Author: Someone
Date: Some day
Place: Some place
Email: mashadab@utexas.edu
"""

#Description: To solve steady unconfined aquifer in Danube-Tiscza Interfluve heterogeneous

#Import functions
#Python
import sys
import numpy as np
import matplotlib.pyplot as plt

#Inhouse
sys.path.insert(1,'../../HW1_BuildGrid_Ops/Python/')  #importing path
sys.path.insert(1,'../../HW2_BuildBnd_SolveLBVP/Python/')
sys.path.insert(1,'../../HW3_CompMean/Python/')

from classfun import *   #For classes and attributes
from build_gridfun2D import build_grid
from build_opsfun2D import build_ops
from build_bndfun_optimized import build_bnd
from solve_lbvpfun_optimized import solve_lbvp
from comp_mean_matrix import comp_mean

#Physical properties
cm2m = 1/100                 #cm to m conversion
yr2s = 365.25 * 24 * 60 * 60 #year to s conversion

Length = 85070               #Distance between Danube and Tiscza rivers [m]
Width  = 5430                #Width of segment considered [m]
K1     = 2e-1*cm2m           #Hydraulic conductivity of the left half segment [m/s]
K2     = 2e-3*cm2m           #Hydraulic conductivity of the right half segment [m/s]
qp     = 1.5*cm2m/yr2s       #Average annual precipitation [m/s]
hD     = 90                  #Height of Danube river [m]
hT     = 80                  #Height of Tiscza river [m]
b      = 100                 #Aquifer thickness [m]

##Generate grid and operators
#Grid
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 100
Grid   = build_grid(Grid)

#Operators
[D,G,C,I,M]  = build_ops(Grid)

#Generate conductivity field
K  = K1*np.ones((Grid.Nx,1)); K[int(Grid.Nx/2):,:] = K2
Kd = comp_mean(K, M, -1, Grid, 1)

L  = - D @ Kd @ G
fs =   qp/(b) * np.ones((Grid.N,1))


##Implement boundary conditions
BC.dof_dir   = np.array([[Grid.dof_xmin],[Grid.dof_xmax]])
BC.dof_f_dir = np.array([[Grid.dof_f_xmin],[Grid.dof_f_xmax]])
BC.g         = np.array([[hD], [hT]])

BC.dof_neu   = np.array([])
BC.dof_f_neu = np.array([])
BC.qb        = np.array([])
#dof_neu not required

[B,N,fn]     = build_bnd(BC,Grid,I) #Building boundary conditions

##Solve linear boundary value problem L * h = fs
u   = solve_lbvp(L, fs+fn, B, BC.g, N)

##Plot solutions
plt.figure()
plt.plot(Grid.xc/1e3,u)
plt.xlabel(r'Horizontal distance, $x$ [km]')
plt.ylabel(r'Water table height, $h$ [m]')
plt.tight_layout()

plt.figure()
plt.semilogy(Grid.xc/1e3,K)
plt.xlabel(r'Horizontal distance, $x$ [km]')
plt.ylabel(r'Hydraulic conductivity, $K$ [m/s]')
plt.tight_layout()

