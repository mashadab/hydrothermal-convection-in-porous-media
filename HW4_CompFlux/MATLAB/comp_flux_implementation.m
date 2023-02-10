%% Define domain and grid
Grid.xmin = 0; 
Grid.xmax = 50; 
Grid.Nx   = 100;
Grid = build_grid(Grid);

%% Define discrete operators
[D,G,C,I,M] = build_ops(Grid);
L  = -D*G+I;                            % system matrix
fs = Grid.xc;                           % r.h.s.
flux = @(h) -G*h;                       % flux
res = @(h,cell) L(cell,:)*h - fs(cell); % residual

%% Define BC's
BC.dof_dir   = [];
BC.dof_f_dir = [];
BC.g         = [];
BC.dof_neu   = [];
BC.dof_f_neu = [];
BC.qb        = [];

[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve boundary value problem

hD = solve_lbvp(L,fs+fn,B,BC.g,N);
qD = comp_flux_gen(flux,res,hD,Grid,BC);


%% Plot results
subplot 211
semilogy(Grid.xc/1e3,hD)
xlabel 'x', ylabel 'h'
%ylim([1e-5 1e-2])
%pbaspect([1 .5 1])

subplot 212
plot(Grid.xf/1e3,qD)
xlabel 'q', ylabel 'x'
%ylim([50 150])
%pbaspect([1 .5 1])