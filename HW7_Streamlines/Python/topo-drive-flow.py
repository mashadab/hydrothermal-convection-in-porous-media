#Topography driven flow
import sys
sys.path.insert(1, '../../HW1_BuildGrid_Ops/Python/')
sys.path.insert(1, '../../HW2_BuildBnd_SolveLBVP/Python/')
sys.path.insert(1, '../../HW3_CompMean/Python/')
sys.path.insert(1, '../../HW4_CompFlux/Python/')
sys.path.insert(1, '../../HW5_Radial_Coordinates/Python/')
sys.path.insert(1, '../../HW6_2D_buildgird_ops/Python/')
sys.path.insert(1, '../../HW7_Streamlines/Python/')



#Import libraries
#Python
import numpy as np
import matplotlib.pyplot as plt

#Inhouse
from classfun import *
from build_gridfun2D import build_grid
from build_opsfun2D_latest import build_ops
from comp_mean_matrix import comp_mean
from build_bndfun_optimized import build_bnd
from solve_lbvpfun_optimized import solve_lbvp
from comp_fluxfun_gen import comp_flux_gen
from comp_streamfun import comp_streamfun

# Physical problem parameters
Length = 200;    # [m] - width of the valley
dh = 15;         # [m] - 1/2 depth of the valley
Height = 50;     # [m] - aquifer thickness
K = 2e-7;        # [m/s] - background conductivity

# Analytic solution
hana = lambda x,z: Height + dh * np.cos(2*np.pi*x/Length)*np.cosh(2*np.pi*z/Length)/np.cosh(2*np.pi*Height/Length);


#Build grid
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 4
Grid.ymin = 0; Grid.ymax = Height; Grid.Ny = 3 
Grid      = build_grid(Grid)


#Build ops
[D,G,C,I,M] = build_ops(Grid)

K           = K*np.ones(Grid.N)  #hydraulic conductivity
Kd          = comp_mean(K, M, -1, Grid, 1)
L           = -D @ Kd @ G 
fs          = np.zeros((Grid.N,1)) #rhs

flux        = lambda h: -Kd @ (G @ h)
res         = lambda h, cell: L[cell] @ h - fs[cell]


#Boundary conditions
BC.dof_dir   = np.array([Grid.dof_ymax])
BC.dof_f_dir = np.array([Grid.dof_f_ymax])
BC.g         = np.transpose([hana(Grid.xc,Grid.yc[-1])])
BC.dof_neu   = np.array([])
BC.dof_f_neu = np.array([])
BC.qb        = np.array([])


[B,N,fn]     = build_bnd(BC, Grid, I)


#Solve LBVP and compute flux
h    = solve_lbvp(L, fs+fn, B, BC.g, N)
q    = comp_flux_gen(flux, res, h, Grid, BC) 


#Calculate stream function
PSI,psi_min,psi_max  = comp_streamfun(q, Grid)

Xc,Yc  = np.meshgrid(Grid.xf,Grid.yf)
plt.contourf(Xc,Yc,PSI,level=20)   



