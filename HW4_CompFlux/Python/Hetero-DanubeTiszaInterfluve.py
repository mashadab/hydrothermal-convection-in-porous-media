import sys
sys.path.insert(1, '../../HW1-Numerics/Python/')
sys.path.insert(1, '../../HW2-BC_LBVP/Python/')

from classfun import *
from build_gridfun2D import build_grid 
from build_opsfun2D import build_ops
from build_bndfun_optimized import build_bnd
from solve_lbvpfun_optimized import solve_lbvp
from comp_mean_matrix import comp_mean
        
#Parameters of the problem
cm2m = 1/100        #cm to metres
yr2s = 365.25*24*60*60 #year to second conversion
Length = 85070      # Distance between Danube and Tisza rivers [m]
Width = 5430        # Width of segment considered [m]
K1 = 2e-1*cm2m      # Hydraulic conductivity [m/s]
K2 = 2e-3*cm2m      # Hydraulic conductivity [m/s]
qp = 1.5*cm2m/yr2s  # Average annual precipitation [m3/m2/s]
hD = 90             # Elevation of Danube river[m]
hT = 80             # Elevation of Tisza river [m]
b = 100             # Aquifer thickness [m]

#building grid and operators
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 100

Grid = build_grid(Grid)

K  = K2*np.ones(Grid.N) #line array
K[Grid.xc<Length/2] = K1
[D,G,C,I,M] = build_ops(Grid)
Kd   = comp_mean(K, M, -1, Grid, 1)
L    = -D @ Kd @ G #Laplacian operator
fs   = qp/b*np.ones((Grid.N,1))

#Boundary conditions
BC.dof_dir   = np.array([Grid.dof_xmin,Grid.dof_xmax])
BC.dof_f_dir = np.array([Grid.dof_f_xmin,Grid.dof_f_xmax]) 
BC.g         = np.array([hD,hT]) 
BC.dof_neu   = np.array([])
BC.dof_f_neu = np.array([]) 
BC.qb        = np.array([])

B,N,fn = build_bnd(BC,Grid,I)

#Solve LBVP
h      = solve_lbvp(L, fs + fn, B, BC.g, N)

#Plot results
plt.figure(figsize=(10,7.5))
plt.plot(Grid.xc/1e3,K/cm2m)
plt.xlabel(r'$x$ [km]')
plt.ylabel(r'$K$ [cm/s]')

plt.figure(figsize=(10,7.5))
plt.plot(Grid.xc/1e3,h)
plt.xlabel(r'$x$ [km]')
plt.ylabel(r'$h$ [m]')