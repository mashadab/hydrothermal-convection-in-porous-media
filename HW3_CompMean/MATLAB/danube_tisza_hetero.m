%% Physical properties
cm2m = 1/100;        % cm to m conversion
yr2s = 365*24*60^2;  % yr to s conversion
Length = 85070;      % Distance between Danube and Tisza rivers [m]
Width = 5430;        % Width of segment considered [m]
K1 =2e-1*cm2m;       % Hydraulic conductivity [m/s]
K2 =2e-3*cm2m;       % Hydraulic conductivity [m/s]
qp = 1.5*cm2m/yr2s;  % Average annual precipitation [m3/m2/s]
hD = 90;             % Elevation of Danube river[m]
hT = 80;             % Elevation of Tisza river [m]
b = 100;             % Aquifer thickness [m]

%% Generate grid and operators
Grid.xmin = 0; 
Grid.xmax = Length; 
Grid.Nx = 100;
Grid = build_grid(Grid);
[D,G,C,I,M] = build_ops(Grid);
% Conductivity field
K = [K1*ones(Grid.Nx/2,1);K2*ones(Grid.Nx/2,1)];
Kd = comp_mean(K,M,-1,Grid,1);
L = -D*Kd*G;
fs = qp/b*ones(Grid.Nx,1);

%% Set boundary conditions
BC.dof_dir = [Grid.dof_xmin;Grid.dof_xmax];
BC.dof_f_dir = [Grid.dof_f_xmin;Grid.dof_f_xmax];
BC.g = [hD;hT];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve problem
h = solve_lbvp(L,fs+fn,B,BC.g,N);

%% Plot results
subplot 211
semilogy(Grid.xc/1e3,K)
xlabel 'x [km]', ylabel 'K [m/s]'
ylim([1e-5 1e-2])
pbaspect([1 .5 1])

subplot 212
plot(Grid.xc/1e3,h)
xlabel 'x [km]', ylabel 'h [m]'
ylim([50 150])
pbaspect([1 .5 1])
