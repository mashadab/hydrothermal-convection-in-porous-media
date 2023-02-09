Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
Grid = build_grid(Grid);
[D,G,C,I,M] = build_ops(Grid);
K = Grid.xc; 
Kd = comp_mean(K,M,1,Grid,1);