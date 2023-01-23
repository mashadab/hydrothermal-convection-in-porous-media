function [D,G,C,I,M]=build_ops(Grid)
% author: 
% date: 
% description:
% This function computes the discrete operator matrices on a
% staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix
% C = discrete curl - not defined in 1D
% I = Nx by Nx identity matrix
% M = Nx+1 by Nx discrete mean matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,C,I,M]=build_ops(Grid);

Nx = Grid.Nx; Nfx = Grid.Nfx;
e = ones(Nx,1);
%% 1) Build sparse Divergence operator
D = 1/Grid.dx*spdiags([-e e],[0 1],Grid.Nx,Grid.Nx+1);

%% 2) Build sparse Gradient operator
%  Interior
G = -D';
% Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax]; % all dof's on boundary
G(dof_f_bnd,:) = 0;

% 3) Discrete Curl operator (not defined in 1D)
C = [];

% 4) Sparse Identity 
I = speye(Nx);

% 5) Sparse Mean
% Interior
M = Grid.dx/2 * abs(G);
% Boundaries
M(1,1) = 1; 
M(Grid.Nx+1,Grid.Nx) = 1;