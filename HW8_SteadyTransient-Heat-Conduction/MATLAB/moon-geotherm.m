'''
Our understanding of the Moon has been revolutionized by NASA's GRAIL Mission. One important insight is that the continuedasteroid impacts onto the surface have created a deep fractured/porous layer with porosities up to 25%. In this exercise we will explore the effect of this near surface porosity on the geotherm of the Moon. We assume the following exponential decrease in porosity with depth
, where  is the surface porosity,  km is the radius of the Moon and  km is the decay depth. Hence you should solve the following heat conduction problem:
 on  with the boundary conditions  and ,
where  kg/m is the mean density,  W/kg is the heat production and  K is the mean surface temperature. The mean thermal conductivity is given by
, where  W/(m K) is the thermal conductivity of the rock. Use the discrete operatore in spherical radial geometry and the appropriately averaged conductivity matrix to solve this problem numerically. Use a numerical gird of Nx = 300.
Note: The comp_mean function used here has 4 arguments! Please set the 3rd argument to 1.
Kd = comp_mean(K,p,1,Grid)

'''


H   = 7.38e-12;% mean heat production [W/kg]
rho = 3.3e3;   % mean density [kg/m^3]
R   = 1738e3;  % Radius of the moon [m]
T0  = 250;     % mean surface temperature [K]
kappa = 3.3;   % mean thermal conductivity [W/(m K)]
phi0  = 0.25;  % surface porosity
hr    = 30e3;  % decay depth [m]

% Porosity decay
phi = @(r) phi0 .* exp(-(R-r)/hr);

% Grid and ops
Grid.xmin = 0; 
Grid.xmax = R; 
Grid.Nx   = 300;
Grid.geom = 'spherical_r';
Grid      = build_grid(Grid);
[D,G,~,I,M] = build_ops(Grid);
kappa_eff = (1-phi(Grid.xc)) * kappa;
size(kappa_eff)
%Kd = comp_mean(kappa_eff, 0, 1, Grid);
Kd = spdiags(M*kappa_eff,0,Grid.Nf,Grid.Nf);
% PDE ops
L = - D * Kd * G; 
fs =  (1 - phi(Grid.xc))* rho * H;
flux = @(u) - Kd * G * u;
res  = @(u,cell) L(cell,:) * u - fs(cell);

BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = T0; 
BC.dof_neu = []; 
BC.dof_f_neu = []; 
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

u = solve_lbvp(L,fs+fn,B,BC.g,N);
q = comp_flux_gen(flux,res,u,Grid,BC);

subplot 411
plot(Grid.xc/1e3,fs*1e6,'b-','markerfacecolor','w','markersize',6)
ylabel '(1-\phi)\rho H [\muW/m3]'

subplot 412
plot(Grid.xc/1e3,kappa_eff,'b-','markerfacecolor','w','markersize',6)
ylabel '(1-\phi)\kappa [W/(m K)]'
ylim([2 4])

subplot 413
plot(Grid.xf/1e3,q*1e3,'b-','markerfacecolor','w','markersize',6), hold off
ylabel 'q [mW/m^2]'
ylim([0 15])
subplot 414
plot(Grid.xc/1e3,u,'b-','markerfacecolor','w','markersize',6)
xlabel 'r [km]', ylabel 'T [K]'
 