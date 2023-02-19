import scipy.sparse as sp
import numpy as np

def build_ops(Grid):
    # author: Mohammad Afzal Shadab
    # date: 1/27/2020
    # description:
    # This function computes the discrete divergence and gradient matrices on a
    # regular staggered grid using central difference approximations. The
    # discrete gradient assumes homogeneous boundary conditions.
    # Input:
    # Grid = structure containing all pertinent information about the grid.
    # Output:
    # D = discrete divergence matrix
    # G = discrete gradient matrix
    # I = identity matrix

    Nx = Grid.Nx
    Ny = Grid.Ny
    N  = Grid.N

    # Two dimensional divergence    
    #     Readable implementation
    #     # 2D divergence matrices
    
    if (Nx>1) and (Ny>1): #2D case
        #One diamentional divergence
        Dy = sp.spdiags(([-np.array(np.ones((Ny+1),'float64')) , np.array(np.ones((Ny+1),'float64'))])/np.asarray(Grid.dy),np.array([0,1]),Ny,Ny+1).tocsr()#.toarray() # Dy^1
        
        #Two dimensional divergence
        Dy = sp.kron(np.eye(Nx), Dy) #y component Dy^2
        
        e  = np.array(np.ones(Ny*(Nx+1),'float64'))
        Dx = sp.spdiags(([-e , e])/np.asarray(Grid.dx),np.array([0,Ny]),N,(Nx+1)*Ny).tocsr()#.toarray() # 2D div-matrix in x-dir

        #D  = np.concatenate((Dx , Dy), axis=1)
        D  = sp.hstack([Dx , Dy])    
        dof_f_bnd = np.concatenate(np.array([Grid.dof_f_xmin-1, Grid.dof_f_xmax-1, Grid.dof_f_ymin-1, Grid.dof_f_ymax-1]), axis=0)       # boundary faces
        dof_f_bnd = np.transpose(dof_f_bnd)
        
    elif (Nx > 1) and (Ny == 1): #one dimensional in x direction
        D = sp.spdiags(([-np.array(np.ones((Nx+1),'float64')),np.array(np.ones((Nx+1),'float64'))])/np.asarray(Grid.dx),np.array([0,1]),Nx,Nx+1).tocsr()#.toarray() # 1D div-matrix in x-dir
        dof_f_bnd = [Grid.dof_f_xmin-1, Grid.dof_f_xmax-1] # boundary faces
        dof_f_bnd = np.transpose(dof_f_bnd)  

    elif (Nx == 1) and (Ny > 1): #one dimensional in y direction
        D = sp.spdiags(([-np.array(np.ones((Ny+1),'float64')),np.array(np.ones((Ny+1),'float64'))])/np.asarray(Grid.dy),np.array([0,1]),Ny,Ny+1).tocsr()#.toarray() # 1D div-matrix in y-dir
        dof_f_bnd = [Grid.dof_f_ymin-1, Grid.dof_f_ymax-1] # boundary faces
        dof_f_bnd = np.transpose(dof_f_bnd)  

    # Gradient
    # Note this is only true in cartesian coordinates!
    # For more general coordinate systems it is worth
    # assembling G and D seperately.
    #print(D)

    D =  sp.csr_matrix(D)
    G = -sp.csr_matrix.transpose(D)
    G =  zero_rows(G,dof_f_bnd)

    #Identity
    I = (sp.eye(Grid.N)).tocsr()
    
    #Curl matrix 
    C = []
    
    #Algebraic mean  
    if Grid.Nx>1 and Grid.Ny==1:
        #Averaging in x-direction considering the zero-flux at boundary
        Avg_x1 = sp.spdiags(([np.array(0.5*np.ones((Grid.Nx),'float64')),0.5*np.array(np.ones((Grid.Nx),'float64'))]),np.array([0,1]),Grid.Nx-1,Grid.Nx)
        Avg_x1 = sp.vstack([np.zeros((1,Grid.Nx)),Avg_x1, np.zeros((1,Grid.Nx))])         
        M      = Avg_x1.copy()
        
    elif Grid.Nx==1 and Grid.Ny>1:
        Avg_y1 = sp.spdiags(([np.array(0.5*np.ones((Grid.Ny),'float64')),0.5*np.array(np.ones((Grid.Ny),'float64'))]),np.array([0,1]),Grid.Ny-1,Grid.Ny)
        Avg_y1 = sp.vstack([np.zeros((1,Grid.Ny)),Avg_y1,np.zeros((1,Grid.Ny))]) 
        M      = Avg_y1.copy()     
    
    elif Grid.Nx>1 and Grid.Ny>1:
        #Averaging in y-direction considering the zero-flux at boundary
        Avg_y1 = sp.spdiags(([np.array(0.5*np.ones((Grid.Ny),'float64')),0.5*np.array(np.ones((Grid.Ny),'float64'))]),np.array([0,1]),Grid.Ny-1,Grid.Ny)
        Avg_y1 = sp.vstack([np.zeros((1,Grid.Ny)),Avg_y1,np.zeros((1,Grid.Ny))]) 
        Avg_y2 = sp.kron(sp.eye(Grid.Nx),Avg_y1)
    
        #Averaging in x-direction considering the zero-flux at boundary
        Avg_x1 = sp.spdiags(([np.array(0.5*np.ones((Grid.Nx),'float64')),0.5*np.array(np.ones((Grid.Nx),'float64'))]),np.array([0,1]),Grid.Nx-1,Grid.Nx)
        Avg_x1 = sp.vstack([np.zeros((1,Grid.Nx)),Avg_x1, np.zeros((1,Grid.Nx))]) 
        Avg_x2 = sp.kron(Avg_x1,sp.eye(Grid.Ny))         
        M      = sp.vstack([Avg_x2, Avg_y2])
 
    if Grid.geom == 'cartesian':
        Grid.A  = np.concatenate([np.ones((Grid.Nfx,1))*Grid.dy,np.ones((Grid.Nfy,1))*Grid.dx,[[Grid.dx*Grid.dy]]], axis=0 )
        Grid.V  = np.ones((Grid.N,1))*Grid.dx*Grid.dy           

    elif Grid.geom == 'cylindrical':
        # In cylindrical coordinates we assume dz gives the height of the cylinder. We don't use dz.
        # assumes x = r
        print('Operators built for cylindrical r-z geometry\n.')

        DOF = np.transpose(np.reshape(Grid.dof,(Grid.Nx,Grid.Ny)))   
        DOFx = np.transpose(np.array([list(range(1,Grid.Nfx+1,1))]).reshape((Grid.Nx+1,Grid.Ny)))
        DOFy = np.transpose(np.reshape(Grid.Nfx + np.array([list(range(1,Grid.Nfy+1,1))]),(Grid.Nx,Grid.Ny+1)))
        R_rfaces = np.transpose([Grid.xf[np.argwhere(DOFx.T)[:,0]]])
        R_zfaces = np.transpose([Grid.xc[np.argwhere(DOFy.T)[:,0]]])
        R_Vol    = np.transpose([Grid.xc[np.argwhere(DOF.T)[:,0]]])
    
        Rcinv = sp.spdiags(1/(R_Vol.T),0,Grid.N,Grid.N)
        
        print(np.shape(Rcinv),np.shape(D[:,:Grid.Nfx]),np.shape(sp.spdiags(R_rfaces.T,0,Grid.Nfx,Grid.Nfx)))
        
        D[:,:Grid.Nfx] =  Rcinv @ D[:,:Grid.Nfx] @ sp.spdiags(R_rfaces.T,0,Grid.Nfx,Grid.Nfx)
        
        #no need to do anything in z direction
        
        
    elif Grid.geom == 'cylindrical_r':
        # In cylindrical coordinates dz are not used.
        # assumes x = r    
        print('Operators built for 1D cylindrical geometry\n.')
        Rf    = sp.spdiags(Grid.xf,0,Grid.Nx+1,Grid.Nx+1)
        Rcinv = sp.spdiags(1/Grid.xc,0,Grid.Nx,Grid.Nx)
        D     = Rcinv @ D @ Rf
    
    elif Grid.geom == 'spherical_r':
        # In spherical coordinates dy and dz are not used.
        # assumes x = r    

        print('Operators built for 1D spherical geometry\n.')
        Rf    = sp.spdiags(Grid.xf**2,0,Grid.Nx+1,Grid.Nx+1)
        Rcinv = sp.spdiags(1/Grid.xc**2,0,Grid.Nx,Grid.Nx)
        D     = Rcinv @ D @ Rf

    else:
        raise ValueError('Unknown grid geometry specified')
 
    return D,G,C,I,M;

def zero_rows(M, rows_to_zero):

    ixs = np.ones(M.shape[0], int)
    ixs[rows_to_zero] = 0
    D = sp.diags(ixs)
    res = D * M
    return res

from build_gridfun_Radial import build_grid
class Grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx   = []

Grid.xmin = 0; Grid.xmax = 4; Grid.Nx = 4
Grid.ymin = 0; Grid.ymax = 3; Grid.Ny = 3
Grid.geom = 'cylindrical'
Grid = build_grid(Grid)
[D,G,C,I,M]  = build_ops(Grid)


