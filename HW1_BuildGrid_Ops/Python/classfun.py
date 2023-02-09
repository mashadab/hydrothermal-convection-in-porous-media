#import libraries
#inbuilt
import numpy as np
import matplotlib.pyplot as plt

#user defined 
import warnings
warnings.filterwarnings('ignore')
from IPython import get_ipython
get_ipython().magic('reset -sf') #for clearing everything
get_ipython().run_line_magic('matplotlib', 'qt') #for plotting in separate window

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'font.family': 'Times'})


#Classes implemented in the code

class Grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx = []

class BC:
    def __init__(self):
        self.dof_dir = []       # identify cells on Dirichlet bnd
        self.dof_f_dir = []     # identify faces on Dirichlet bnd
        self.dof_neu = []       # identify cells on Neumann bnd
        self.dof_f_neu = []     # identify faces on Neumann bnd
        self.g = []             # column vector of non-homogeneous Dirichlet BCs (Nc X 1)
        self.qb = []         
        
    class h:
        def __init__(self):
            self.dof_dir = []       # identify cells on Dirichlet bnd
            self.dof_f_dir = []     # identify faces on Dirichlet bnd
            self.dof_neu = []       # identify cells on Neumann bnd
            self.dof_f_neu = []     # identify faces on Neumann bnd
            self.g = []             # column vector of non-homogeneous Dirichlet BCs (Nc X 1)
            self.qb = []     
            
    class u:
        def __init__(self):
            self.dof_dir = []       # identify cells on Dirichlet bnd
            self.dof_f_dir = []     # identify faces on Dirichlet bnd
            self.dof_neu = []       # identify cells on Neumann bnd
            self.dof_f_neu = []     # identify faces on Neumann bnd
            self.g = []             # column vector of non-homogeneous Dirichlet BCs (Nc X 1)
            self.qb = []     
        
    class phi:
        def __init__(self):
            self.dof_dir = []       # identify cells on Dirichlet bnd
            self.dof_f_dir = []     # identify faces on Dirichlet bnd
            self.dof_neu = []       # identify cells on Neumann bnd
            self.dof_f_neu = []     # identify faces on Neumann bnd
            self.g = []             # column vector of non-homogeneous Dirichlet BCs (Nc X 1)
            self.qb = []     
                
        
        
        