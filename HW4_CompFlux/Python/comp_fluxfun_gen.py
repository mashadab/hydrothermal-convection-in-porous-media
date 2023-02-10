# import python libraries
import numpy as np

def comp_flux_gen(flux,res,u,Grid,Param):
      # author: Mohammad Afzal Shadab
      # date: 11 February 2021
      # Description:
      # Computes the fuxes on the interior from flux(u) and reconstructs the
      # fluxes on the boundary faces from the residuals in the adjacent boundary
      # cells, res(u). 
      #
      # Note to self: This is an attempt to make the function more general
      #               the real test will be if this works with variable
      #               densities.
    
      # Input:
      # flux = anonymous function computing the flux (correct in the interior)
      # res = anonymous function computing the residual 
      # u = vector of 'flux potential' (head, temperature,electric field,...)
      # Grid = structure containing pertinent information about the grid
      # Param = structure containing pertinent information about BC's
      # 
      # Output:
      # q = correct flux everywhere
    
    ## Compute interior fluxes
    q = flux(u)
    
    ## Compute boundary fluxes
    #1) Identify the faces and cells on the boundary
    # note: check if o.k. for homogeneous Neumann problem
    if not Param.dof_neu.any(): 
        if Grid.Ny == 1:
            dof_cell = Param.dof_dir
            dof_face = Param.dof_f_dir
        else:       
            dof_cell = (np.squeeze(Param.dof_dir, axis=0))
            dof_face = (np.squeeze(Param.dof_f_dir, axis=0))

    elif not Param.dof_dir.any(): 
        dof_cell = Param.dof_neu
        dof_face = Param.dof_f_neu
    
    #For non-empty Dirichlet and Neumann BC
    else:
        if Grid.Ny > 1:
            dof_cell = np.concatenate((np.squeeze(Param.dof_dir, axis=0),np.squeeze(Param.dof_neu, axis=0)),axis=0)
            dof_face = np.concatenate((np.squeeze(Param.dof_f_dir, axis=0),np.squeeze(Param.dof_f_neu, axis=0)),axis=0)
        else:
            dof_cell = np.concatenate((Param.dof_dir,Param.dof_neu),axis=0)
            dof_face = np.concatenate((Param.dof_f_dir,Param.dof_f_neu),axis=0)
        
    dof_cell = list(filter(None, dof_cell))
    dof_face = list(filter(None, dof_face))
    
    # 2) Determine sign of flux: Convention is that flux is positive in
    #    coordinate direction. So the boundary flux, qb is not equal to q*n,
    #    were n is the outward normal!
    sign = np.multiply(np.isin(dof_face,np.concatenate((Grid.dof_f_xmin,Grid.dof_f_ymin),axis=0)), 1) - \
           np.multiply(np.isin(dof_face,np.concatenate((Grid.dof_f_xmax,Grid.dof_f_ymax),axis=0)), 1)
     
    #Because of Python indexing
    dof_cell = np.subtract(dof_cell,1)
    dof_face = np.subtract(dof_face,1)
    
    # 3) Compute residuals and convert them to bnd fluxes    
    q[dof_face,:] =  np.transpose([sign]) * res(u,dof_cell) *Grid.V[dof_cell,:]/Grid.A[dof_face,:]

    return q;