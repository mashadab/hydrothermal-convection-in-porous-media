
%{
Here we compute the flow net for topography-driven flow beneath a simple sinusoudal topography.
This is a storied topic going back to a seminal paper by Toth in 1940 and you can read more about
it in this free text-book. 
The ideas is the following. In areas with sufficient precipitation and moderat topography, the water
table is a muted image of the topography, as shown in the image from Hubbert (1940), below.

Then it is generally assumed that the depth of the aquifer is large with respect to the amplitude of
of the topographic variation, so that the simulation domain remains a rectangle. Note that this assumption
does not hold fo the image above.

The simplified mathematical problem is given by:
PDE:  on  and 
BC's:  and  elsewhere
This problem has the analytic solution:
.
Use the analytic solution to set the boundary condition in the cell centers along the top boundary.
%}

%% Physical problem parameters
Length = 200;     % [m] - width of the valley
dh = 15;          % [m] - 1/2 depth of the valley
Height = 50;      % [m] - aquifer thickness
K = 2e-7;        % [m/s] - background conductivity

%% Analytic solution
hana = @(x,z) Height + dh * cos(2*pi*x/Length).*cosh(2*pi*z/Length)/cosh(2*pi*Height/Length);

%% Set up grid
Grid.xmin = 0; 
Grid.xmax = Length; 
Grid.Nx = 80;
Grid.ymin = 0; 
Grid.ymax = Height; 
Grid.Ny = 20;

Grid    = build_grid2D(Grid); 

%% Set boundary conditions
BC.dof_dir   = Grid.dof_ymax;
BC.dof_f_dir = Grid.dof_f_ymax;
BC.g         = hana(Grid.xc,Grid.yc(end));
BC.dof_neu   = [];
BC.dof_f_neu = [];
BC.qb        = [];

%% Set permeability field
K = K*ones(Grid.N,1); % column vector of cell center hydraulic conductivities

%% Discrete differential operators
[D,G,C,I,M]  = build_ops2D(Grid);
Kd           = comp_mean(K,M,1,Grid,1);
L            = - D * Kd * G;
fs           = spalloc(Grid.N,1,0);

flux         = @(h) -Kd * G * h;  
res          = @(h,cell) L(cell,:) * h - fs(cell);

%% Build boundary operators
[B,N,fn]     = build_bnd(BC,Grid,I);

%% Solve system
h = solve_lbvp(L,fs+fn,B,BC.g,N); 
q = comp_flux_gen(flux,res,h,Grid,BC);

%% Run learner solution.
[PSI, psi_min, psi_max] = comp_streamfun(q,Grid);

subplot 211
plot(Grid.xf,hana(Grid.xf,Height))
pbaspect([1 .2 1])
xlabel 'x [m]'
ylabel 'h_b(x) [m]'

subplot 212
plot_flownet(10,10,h,PSI,'b-','r-',Grid)
