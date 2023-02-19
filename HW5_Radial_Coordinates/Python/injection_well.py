
##Import libraries
#Python
import numpy as np

#Inhouse
import sys
sys.path.insert(1,'../../HW1_BuildGrid_Ops/Python')
sys.path.insert(1,'../../HW2_BuildBnd_SolveLBVP/Python')
sys.path.insert(1,'../../HW3_CompMean/Python')
sys.path.insert(1,'../../HW4_CompFlux/Python')
sys.path.insert(1,'../../HW5_Radial_Coordinates/Python')

from classfun import *
from build_gridfun_Radial import build_grid
from build_opsfun_Radial import build_ops
from build_bndfun_optimized import build_bnd
from solve_lbvpfun_optimized import solve_lbvp
from comp_fluxfun_gen import comp_flux_gen

##Injection well simulation parameters
rw = .1   # well radius
r0 = 100  # domain radius
h0 = 1    # far-field head
Qw = 1    # well flow rate
H = 1     # aquifer thickness
K = 1     # hydraulic conductivity
Aw = 2*np.pi*rw*H # wellbore area
Nx = 10

##Analytic solution
xa = np.linspace(rw,r0,1000)
ha = lambda r: h0 - Qw/(2*np.pi*K*H) * np.log(r/r0) #Head [m]
qa = lambda r: Qw/(2*np.pi*H*r)                     #Vol flux [m/s]
Qa = lambda r: Qw + 0*r                             #Flow rate [m^3/s]

##Numerical Solution

#Fluid flux at the well
qw = Qa(rw)/Aw

#Grid and operator
Grid.xmin = rw; Grid.xmax = r0; Grid.Nx = Nx 
Grid.geom = 'cylindrical_r'


Grid      = build_grid(Grid)
[D,G,C,I,M]= build_ops(Grid)
L         = - K * D @ G
fs        = np.zeros((Grid.N,1))
flux      = lambda h: - K* G @ h
res       = lambda h,dof_cell: L[dof_cell,:] @ h - fs[dof_cell]

#Boundary conditions
BC.dof_dir  = np.array([Grid.dof_xmax])
BC.dof_f_dir= np.array([Grid.dof_f_xmax])
BC.g        = np.array([ha(Grid.xc[Grid.dof_xmax-1])])

BC.dof_neu= np.array([Grid.dof_xmin])
BC.dof_f_neu=np.array([Grid.dof_f_xmin])
BC.qb     = np.array([qw]) 

[B,N,fn]  = build_bnd(BC,Grid,I)

#Compute solution
h         = solve_lbvp(L, fs+fn, B, BC.g, N)
q         = comp_flux_gen(flux, res, h, Grid, BC)
Q         = q*Aw/rw*np.transpose([Grid.xf])

plt.figure(figsize=(15,5))
plt.subplot(131)
plt.plot(xa,ha(xa),'b-',label='Analytical')
plt.plot(Grid.xc,h,'ro',label='Numerical')
plt.legend()
plt.xlabel('r'); plt.ylabel('h')

plt.subplot(132)
plt.plot(xa,qa(xa),'b-')
plt.plot(Grid.xf,q,'ro')
plt.xlabel('r'); plt.ylabel('q')

plt.subplot(133)
plt.ylim([0.5,1.5])
plt.plot(xa,Qa(xa),'b-')
plt.plot(Grid.xf,Q,'ro')
plt.xlabel('r'); plt.ylabel('Q')
plt.tight_layout()

