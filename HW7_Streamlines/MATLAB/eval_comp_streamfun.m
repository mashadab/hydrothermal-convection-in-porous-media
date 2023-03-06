%% Set up grid and operators
Grid.xmin = 0; Grid.xmax = 2; Grid.Nx = 41;
Grid.ymin = 0; Grid.ymax = 2; Grid.Ny = 50;
Grid.psi_x0 = 'xmin_ymin';
Grid.psi_dir = 'yx';
Grid.geom    = 'cartesian';
Grid = build_grid2D(Grid);

%% Set boundary conditions
BC.dof_dir = [Grid.dof_xmin;Grid.dof_ymax(2:Grid.Nx)]; 
BC.dof_f_dir = [Grid.dof_f_xmin;Grid.dof_f_ymax(2:Grid.Nx)];
BC.g = [ones(Grid.Ny,1);zeros(Grid.Nx-1,1)];
BC.dof_neu = []; BC.dof_f_neu = []; BC.qb = [];

%% Set permeability field
K = ones(Grid.N,1);

%% Discrete differential operators
[D,G,~,I,M] = build_ops2D(Grid);
Kd = comp_mean(K,M,-1,Grid,1);
L = -D*Kd*G;  fs = spalloc(Grid.N,1,0);
flux = @(h) -Kd*G*h;
res = @(h,cell) L(cell,:)*h-fs(cell);

%% Build boundary operators
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve system
h = solve_lbvp(L,fs+fn,B,BC.g,N); 
q = comp_flux_gen(flux,res,h,Grid,BC);

%% Compute stream function and plot flownet
[PSI, psi_min, psi_max] = comp_streamfun(q,Grid);
%plot_flownet(11,11,h,PSI,'b-','r',Grid)