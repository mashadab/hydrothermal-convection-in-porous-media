function [D,G,C,I,M]=build_ops_Radial(Grid)
% author: Marc Hesse
% date: 09/08/2014
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix
% C = discrete curl - not defined in 1D
% I = Nx by Nx identity matrix
% M = Nx+1 by Nx discrete mean matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; Grid.geom =
% 'cylindrical_r'
% >> Grid = build_grid_Radial(Grid);
% >> [D,G,C,I,M]=build_ops_Radial(Grid);

Nx = Grid.Nx; Nfx = Grid.Nfx; 
% 1) Build sparse Divergence operator
D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
% 2) Obtain sparse Gradient operator in interior
G = -D';
% ) Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [Grid.dof_f_xmin,Grid.dof_f_xmax]; % all dof's on boundary
G(dof_f_bnd,:) = 0;

% Modify D for geometry
if strcmp(Grid.geom,'cylindrical_r')
    fprintf('Operators built for 1D cylindrical geometry\n.')
    Rf = spdiags(Grid.xf,0,Grid.Nx+1,Grid.Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Grid.Nx,Grid.Nx);
    D = Rcinv * D * Rf;
elseif strcmp(Grid.geom,'spherical_r')
    fprintf('Operators built for 1D spherical geometry\n.')
    Rf = spdiags(Grid.xf.^2,0,Grid.Nx+1,Grid.Nx+1);
    Rcinv = spdiags(1./Grid.xc.^2,0,Grid.Nx,Grid.Nx);
    D = Rcinv * D * Rf;
elseif strcmp(Grid.geom,'cartesian')
    % nothing needs to be done
else
    error('Unknown geometry.')
end

% 3) Discrete Curl operator (not defined in 1D)
C = [];

% 4) Identity
I = speye(Grid.Nx);

% 5) Sparse Mean
% Interior
M = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx);
% Boundaries
M(1,1) = 1; M(Nfx,Nx) = 1;
    