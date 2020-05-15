% Generate the Discrete Laplacians based on French grid

load ../datasets/grid.mat; % load the grid generated with ../datasets/generate_grid.m

J = J2*J1; % total grid point numbers

%% Laplacians for ADI schemes

%% interior 
% Ly is the discretized Laplacian taking indices columnwise (along Y):        
%
%  1   J1+1  2*J1+1  .  .  .  (J2-1)*J1+1
%  2   J1+2  2*J1+2           (J2-1)*J1+2
%  3   J1+3  2*J1+3           (J2-1)*J1+3
%  .   .     .       .        .
%  .   .     .          .     .  
%  .   .     .             .  .  
%  J1  2*J1  .                J2*J1
%


% 1D discretized Laplcien for the column dominant indexes
Ly = sparse(interior,interior,-2,J,J); 
Ly = Ly + sparse(interior,interior+1,1,J,J);
Ly = Ly + sparse(interior,interior-1,1,J,J);

% Building the row dominant indexes
[index_x,index_y]=ind2sub([J1,J2],interior);
interiorX=sub2ind([J2,J1],index_y,index_x);

% 1D discretized Laplcien for the row dominant indexes
Lx = sparse(interiorX,interiorX,-2,J,J);
Lx = Lx + sparse(interiorX,interiorX+1,1,J,J);
Lx = Lx + sparse(interiorX,interiorX-1,1,J,J);

save L2D_ADI Lx Ly interiorX % save the Laplacian in a mat-file L2D_ADI.mat
