# Stokes solver
# author: Mohammad Afzal Shadab
# email: mashadab@utexas.edu
# date: 04/09/2022
import sys
sys.path.insert(1, '../../HW1_BuildGrid_Ops/Python/')
sys.path.insert(1, '../../HW2_BuildBnd_SolveLBVP/Python/')
sys.path.insert(1, '../../HW3_CompMean/Python/')
sys.path.insert(1, '../../HW4_CompFlux/Python/')
sys.path.insert(1, '../../HW5_Radial_Coordinates/Python/')
sys.path.insert(1, '../../HW6_2D_buildgird_ops/Python/')
sys.path.insert(1, '../../HW7_Streamlines/Python/')

# import Python libraries
from scipy.sparse import bmat,csr_matrix
import matplotlib.pyplot as plt

# import personal libraries
from classfun import *  #importing the classes and relevant functions
from build_stokes_grid_fun import build_stokes_grid
from build_stokes_ops_fun import build_stokes_ops
from build_bndfun_optimized import build_bnd
from solve_lbvpfun_optimized import solve_lbvp
from quiver_plot import quiver_plot
from comp_streamfun import comp_streamfun
from plot_streamlines import plot_streamlines

#problem parameter
mu  = 1.0  #Viscosity nondimensionalized 

#building grid
Gridp.xmin = 0.0 ; Gridp.xmax = 1 ; Gridp.Nx   = 100
Gridp.ymin = 0.0 ; Gridp.ymax = 1 ; Gridp.Ny   = 100
Grid = build_stokes_grid(Gridp)

#simulation name
simulation_type = 'lid_driven_cavity_flow_with_no_slip'   #lid_driven_cavity_flow_with_slip or 'no_flow' 
simulation_name = f'stokes_solver_test{simulation_type}_domain{Gridp.xmax-Gridp.xmin}by{Gridp.ymax-Gridp.ymin}_N{Gridp.Nx}by{Gridp.Ny}'

#building operators
D, Edot, Dp, Gp, I = build_stokes_ops(Grid)

A  = 2.0 * mu * D @ Edot
L  = bmat([[A, -Gp], [Dp, None]],format="csr")
fs = csr_matrix((Grid.N, 1), dtype=np.float64)



#Boundary conditions
if 'lid_driven_cavity_flow_with_slip' in simulation_type:
    BC.dof_dir =  np.concatenate((Grid.dof_pene, \
                                  Grid.dof_ymax_vt,\
                                  Grid.dof_pc))
    BC.g = np.transpose([np.concatenate((np.zeros(Grid.N_pene), \
                                         np.ones(len(Grid.dof_ymax_vt)),\
                                        [0.0]))])
    
elif 'lid_driven_cavity_flow_with_no_slip' in simulation_type: #did not work
    BC.dof_dir =  np.concatenate((Grid.dof_pene, \
                                  Grid.dof_ymax_vt[1:-1],\
                                  Grid.dof_ymin_vt[1:-1],\
                                  Grid.dof_xmax_vt[1:-1],\
                                  Grid.dof_xmin_vt[1:-1],\
                                  Grid.dof_pc))
    BC.g = np.transpose([np.concatenate((np.zeros(Grid.N_pene), \
                                         np.ones(len(Grid.dof_ymax_vt[1:-1])),\
                                         np.zeros(len(Grid.dof_ymin_vt[1:-1])),\
                                         np.zeros(len(Grid.dof_xmax_vt[1:-1])),\
                                         np.zeros(len(Grid.dof_xmin_vt[1:-1])),\
                                        [0.0]))])

else:#no flow
    BC.dof_dir =  np.concatenate((Grid.dof_solid_bnd, \
                                  Grid.dof_pc))
    BC.g        = np.concatenate((np.zeros(Grid.N_solid_bnd), \
                                  [0.0]))  

BC.dof_neu = np.array([])
[B,N,fn] = build_bnd(BC,Grid,I)

#Solving for Stokes flow
u = solve_lbvp(L,fs+fn,B,BC.g,N)
v = u[:Grid.p.Nf,:]; p = u[Grid.p.Nf+1:,:] #Extracting velocity and pressure inside
[PSI,psi_min,psi_max] = comp_streamfun(v,Gridp)

#Plotting
#quiver_plot(simulation_name,Grid,v)
plot_streamlines(simulation_name,Grid,v,PSI,psi_min,psi_max,'label_no')
