function [Grid] = build_grid(Grid)
% Author: Marc Hesse
% Date: 09/12/2014
% author: your name
% date: today
% Description:
% This function computes takes in minimal definition of the computational
% domain and grid and computes all containing all pertinent information 
% about the grid. 
% Input:
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx   = number of grid cells
% 
% Output: (suggestions)
% Grid.Lx = scalar length of the domain
% Grid.dx = scalar cell width
% Grid.Nfx = number of fluxes in x-direction
% Grid.xc = Nx by 1 column vector of cell center locations
% Grid.xf = Nfx by 1 column vector of cell face locations
% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell numbers
% Grid.dof_xmin  = scalar cell degree of freedom corrsponding to the left boundary
% Grid.dof_xmax  = scalar cell degree of freedom corrsponding to the right boundary
% Grid.dof_f_xmin = scalar face degree of freedom corrsponding to the left boundary
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right boundary
% Grid.A = Nfx by 1 vector of the areas of the faces
% Grid.V = Nx by 1 vector of the volues of the cells
% + anything else you might find useful
%
% Example call: 
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% >> Grid = build_grid(Grid);

%% Set up catesian geometry
if ~isfield(Grid,'xmin'); Grid.xmin = 0;  fprintf('Grid.xmin is not defined and has been set to zero.\n');end
if ~isfield(Grid,'xmax'); Grid.xmax = 10; fprintf('Grid.xmax is not defined and has been set to 10.\n'); end
if ~isfield(Grid,'Nx');   Grid.Nx   = 10; fprintf('Grid.Nx is not defined and has been set to 10.\n');end
Grid.Lx = Grid.xmax-Grid.xmin;    % domain length in x
Grid.dx = Grid.Lx/Grid.Nx;        % dx of the gridblocks

if ~isfield(Grid,'ymin'); Grid.ymin = 0; end
if ~isfield(Grid,'ymax'); Grid.ymax = 1; end
if ~isfield(Grid,'Ny');   Grid.Ny   = 1; end
Grid.Ly = Grid.ymax-Grid.ymin;    % domain length in y
Grid.dy = Grid.Ly/Grid.Ny;        % dy of the gridblocks

if ~isfield(Grid,'zmin'); Grid.zmin = 0; end
if ~isfield(Grid,'zmax'); Grid.zmax = 1; end
if ~isfield(Grid,'Nz');   Grid.Nz   = 1; end
Grid.Lz = Grid.zmax-Grid.zmin;    % domain length in z
Grid.dz = Grid.Lz/Grid.Nz;        % dz of the gridblocks

%% Number for fluxes
Grid.Nfx = (Grid.Nx+1);

% Set up mesh
% cell centers 'xc' and cell faces 'xf'   
Grid.xc = [Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax-Grid.dx/2]'; % x-coords of gridblock centers
Grid.xf = [Grid.xmin:Grid.dx:Grid.xmax]'; % x-coords of gridblock faces


%% Set up dof vectors
Grid.dof = [1:Grid.Nx]';              % cell centered degree of freedom/gridblock number
Grid.dof_f = [1:Grid.Nfx]';            % face degree of freedom/face number

%% Boundary dof's
% Boundary cells
Grid.dof_xmin = Grid.dof(1);
Grid.dof_xmax = Grid.dof(Grid.Nx);
% Boundary faces
Grid.dof_f_xmin = Grid.dof_f(1);
Grid.dof_f_xmax = Grid.dof_f(Grid.Nfx);

%% Face areas and cell volumes
Grid.A  = Grid.dy * Grid.dz * ones(Grid.Nx+1,1);
Grid.V  = Grid.dx * Grid.dy * Grid.dz * ones(Grid.Nx,1);