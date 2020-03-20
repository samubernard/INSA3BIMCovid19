%% Generate the Discrete Laplacians based on French grid

load ../datasets/grid.mat; % load the grid generated with ../datasets/generate_grid.m

J = J2*J1; % total grid point numbers

%% interior 
% L is the discretized Laplacian taking indices columnwise (along Y):        
%
%  1   J1+1  2*J1+1  .  .  .  (J2-1)*J1+1
%  2   J1+2  2*J1+2           (J2-1)*J1+2
%  3   J1+3  2*J1+3           (J2-1)*J1+3
%  .   .     .       .        .
%  .   .     .          .     .  
%  .   .     .             .  .  
%  J1  2*J1  .                J2*J1
%

% 2D Discrete Laplacian 
L = sparse(interior, interior,-4,J,J); 
L = L + sparse(interior,interior+1,1,J,J);
L = L + sparse(interior,interior-1,1,J,J);
L = L + sparse(interior,interior+J1,1,J,J);
L = L + sparse(interior,interior-J1,1,J,J);

save L2D L  % save the Laplacian in a mat-file L2D.mat

% see what the discete Laplacian looks like
spy(L)

