%% Generate 2D grid with metropolitan France border

% import image as a 3D-array
A = imread('France_population_density.png');

figure(1); clf;
image(A);
axis equal; % make pixels square

[x,y]=ginput;

save border.mat x y

% set grid size 
A = imread('France_population_density.png');
[ay,ax,~] = size(A);

J1 = 101; % 101 along vertical axis
J2 = round(ax/ay*J1); % keep x:y ratio 
NG = J2*J1; % nbr of points in the grid

% build the grid in matrix order (columnwise) 
%
%  1   J1+1  2*J1+1  .  .  .  (J2-1)*J1+1
%  2   J1+2  2*J1+2           (J2-1)*J1+2
%  3   J1+3  2*J1+3           (J2-1)*J1+3
%  .   .     .       .        .
%  .   .     .          .     .  
%  .   .     .             .  .  
%  J1  2*J1  .                J2*J1
%

%% Find the border 
% Rescale the border (xs,ys) 

xs = x/ax*J2;
ys = y/ay*J1;

plot(xs,ys);
axis ij;
axis equal;

% Define the border
xsint = interp(xs,100);
ysint = interp(ys,100);
border = ceil(ysint) + J1*ceil(xsint); % ceil: make sure the integers are >= 1
border = unique(border);

% Check that the border is completly closed:
G = sparse(J1,J2);
G(border) = 1;
spy(G);

% define the interior: more tricky    
% we will use an iterative method, starting with a seed in the interior
i0 = 50; j0 = 60; % a point in the interior
interior = i0 + J1*j0;
downright = 1:(J2*J1);
upleft = (J2*J1):-1:1;
downleft = repmat((1:J1)',1,J2) + repmat((J2-1)*J1:-J1:0,J1,1); 
upright  = repmat((J1:-1:1)',1,J2) + repmat(0:J1:(J2-1)*J1,J1,1); 

for n = 1:2 % nbr of iterations depends on the shape of the border 
  % a point is in the interior if it is a direct neighbour (up,down,left,right) of a point in the interior
  % and is not a border
  % check all 4 neighbours

  for k = downright
    % check if k has a neighbor in th interior
    if any(interior == (k+1)) && ~any(border == k) % k is above an interior point
      interior = [interior; k];
    end
    if any(interior == (k-1)) && ~any(border == k)  % k is below an interior point
      interior = [interior; k];
    end
    if any(interior == (k+J1))  && ~any(border == k) % k is to the right of an interior point
      interior = [interior; k];
    end
    if any(interior == (k-J1)) && ~any(border == k)  % k is to the left of an interior point
      interior = [interior; k];
    end
  end
  interior = unique(interior);

  for k = upleft
    % check if k has a neighbor in th interior
    if any(interior == (k+1)) && ~any(border == k) % k is above an interior point
      interior = [interior; k];
    end
    if any(interior == (k-1)) && ~any(border == k)  % k is below an interior point
      interior = [interior; k];
    end
    if any(interior == (k+J1))  && ~any(border == k) % k is to the right of an interior point
      interior = [interior; k];
    end
    if any(interior == (k-J1)) && ~any(border == k)  % k is to the left of an interior point
      interior = [interior; k];
    end
  end
  interior = unique(interior);

  
  for k = downleft(:)'
    % check if k has a neighbor in th interior
    if any(interior == (k+1)) && ~any(border == k) % k is above an interior point
      interior = [interior; k];
    end
    if any(interior == (k-1)) && ~any(border == k)  % k is below an interior point
      interior = [interior; k];
    end
    if any(interior == (k+J1))  && ~any(border == k) % k is to the right of an interior point
      interior = [interior; k];
    end
    if any(interior == (k-J1)) && ~any(border == k)  % k is to the left of an interior point
      interior = [interior; k];
    end
  end
  interior = unique(interior);

  for k = upright(:)'
    % check if k has a neighbor in th interior
    if any(interior == (k+1)) && ~any(border == k) % k is above an interior point
      interior = [interior; k];
    end
    if any(interior == (k-1)) && ~any(border == k)  % k is below an interior point
      interior = [interior; k];
    end
    if any(interior == (k+J1))  && ~any(border == k) % k is to the right of an interior point
      interior = [interior; k];
    end
    if any(interior == (k-J1)) && ~any(border == k)  % k is to the left of an interior point
      interior = [interior; k];
    end
  end
  interior = unique(interior);
  
end

G = sparse(J1,J2);
G(interior) = 1;
G(border) = 2;
imagesc(G)

save grid.mat border interior J1 J2






