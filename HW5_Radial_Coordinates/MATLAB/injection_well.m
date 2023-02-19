rw = .1;   % well radius
r0 = 100;  % domain radius
h0 = 1;    % far-field head
Qw = 1;    % well flow rate
H = 1;     % aquifer thickness
K = 1;     % hydraulic conductivity
Aw = 2*pi*rw*H; % wellbore area
Nx = 10;

%% Analytic solution
xa = linspace(rw,r0,1e3);
ha = @(r) h0 - Qw/2/pi/H/K * log(r/r0);
qa = @(r) Qw/2/pi/H./r;
Qa = @(r) Qw+0*r;

%% Numerical solution

% Fluid flux at well
qw = Qa(rw)./Aw;

%% Grid and operators
Grid.xmin = rw; 
Grid.xmax = r0; 
Grid.Nx   = Nx; 
Grid.geom = 'cylindrical_r';
Grid      = build_grid_Radial(Grid);

[D,G,C,I,M] = build_ops_Radial(Grid);
L  = -D*K*G; 
fs = zeros(Grid.Nx,1);
flux = @(h) - K*G*h;
res  = @(h,dof_cell) L(dof_cell,:) * h - fs(dof_cell) ;

%% Boundary conditions
BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g         = ha(Grid.xc(Grid.dof_xmax));

BC.dof_neu   = Grid.dof_xmin;
BC.dof_f_neu = Grid.dof_f_xmin;
BC.qb        = qw;

[B,N,fn] = build_bnd(BC,Grid,I);

%% Compute solution
h = solve_lbvp(L,fs+fn,B,BC.g,N);
q = comp_flux_gen(flux,res,h,Grid,BC);
Q = q.*Aw/rw.*Grid.xf;

subplot 131
plot(xa,ha(xa),'r-'), hold on
plot(Grid.xc,h,'bo')
pbaspect([1 .8 1])
xlabel 'r'
ylabel 'h'
legend('analytic','numeric')

subplot 132
plot(xa,qa(xa),'r-'), hold on
plot(Grid.xf,q,'bo')
pbaspect([1 .8 1])
xlabel 'r'
ylabel 'q'

subplot 133
plot(xa,Qa(xa),'r-'), hold on
plot(Grid.xf,Q,'bo')
pbaspect([1 .8 1])
xlabel 'r'
ylabel 'Q'
ylim([0 1.2])
