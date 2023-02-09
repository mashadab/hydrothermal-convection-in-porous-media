% file: DanubeTiszaInterfluve.m
% date: 8 June 1867
% author: Empress Elisabeth of Austria

% Descrition: Solve for the steady confined aquifer in the Danube Tisza Interfluve

%% Physical properties
cm2m = 1/100;        % cm to m conversion
yr2s = 365*24*60^2;  % yr to s conversion

Length = 85070;      % Distance between Danube and Tisza rivers [m]
Width = 5430;        % Width of segment considered [m]
K =2e-2*cm2m;        % Hydraulic conductivity [m/s]
qp = 1.5*cm2m/yr2s;  % Average annual precipitation [m3/m2/s]
hD = 90;             % Elevation of Danube river[m]
hT = 80;             % Elevation of Tisza river [m]
b = 100;             % Aquifer thickness [m]

%% Generate grid and operators
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 100;
Grid = build_grid(Grid);
[D,G,C,I,M] = build_ops(Grid);
L = -D * G;
fs = qp/(b*K) * ones(Grid.N,1);

%% Set boundary conditions
BC.dof_dir = [Grid.dof_xmin, Grid.dof_xmax]';
BC.g = [hD, hT]';
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve problem
h = solve_lbvp(L,fs+fn,B,BC.g,N);

%% Plot results
plot(Grid.xc/1e3,h)
xlabel 'x [km]', ylabel 'h [m]'
pbaspect([1 .8 1])