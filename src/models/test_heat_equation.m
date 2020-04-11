% Heat equations finite difference implicit Crank-Nicolson
% Domain of integration: Metropolitan France (datasets/grid.mat)
% Discretisation: laplacian/L2D.math
%
% du/dt = D (d^2 u/dx^2 + d^2 u/dx^2)
%
% With Neumann no flux boundary conditions
% Initial conditions: 4 local peaks

% Equation parameters 
D = 400; % diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))

% Spatial grid in grid.mat
% Discretisation of the Laplacian in L2D.mat

load ../datasets/grid.mat
load ../laplacian/L2D.mat
load ../datasets/country.mat

% Simulation parameters 
% The grid has J2 points (in x) by J1 points (in y)
% The French map:
% width: 874 px and 175 px = 200km -> 874 px = 998 km 
% h = 998 km / (J1-1)
% Space units: km 
S  = 998.0;                 % space scale factor;
h = S/(J1-1);               % grid square size
x = h*(0:J2-1);             % space discretisation in x
y = h*(0:J1-1);             % space discretisation in y
[X,Y] = meshgrid(x,y);      % matrix version of (x,y) for plotting purpose
J  = J1*J2;                 % grid size 

% Dynamical variables 
u = zeros(J,1); % u is the current state variable, do not keep memory of past times 
newu = zeros(J,1); 

% Initial Conditions
u( (X(:) - 550).^2 + (Y(:) - 250).^2 < 300 ) = 1; % piecewise constant initial condition
u( (X(:) - 850).^2 + (Y(:) - 270).^2 < 300 ) = 1; % piecewise constant initial condition
u( (X(:) - 670).^2 + (Y(:) - 650).^2 < 300 ) = 1; % piecewise constant initial condition
u( (X(:) - 330).^2 + (Y(:) - 700).^2 < 300 ) = 1; % piecewise constant initial condition

exterior = setdiff(1:J, union(border,interior));

u(exterior) = nan;

X(exterior) = nan;
Y(exterior) = nan;

% time parameters
t0 = 0;
tfinal = 30; 
t = t0;
dt = 0.1;

% camera zoom and orbit
nbr_t = ceil((tfinal-t0)/dt)+1;
dtheta = linspace(-37.5,37.5,floor(nbr_t/2));
dphi = linspace(-10,60,floor(nbr_t/2));
zo = linspace(5,1,floor(nbr_t/2));
dy = linspace(-3,0,floor(nbr_t/2));
dtheta = fliplr([dtheta, repmat(37.5,1,nbr_t - floor(nbr_t/2))]);
dphi = fliplr([dphi, repmat(60,1,nbr_t - floor(nbr_t/2))]);
zo = fliplr([zo, ones(1,nbr_t - floor(nbr_t/2))]);
dy = fliplr([dy, zeros(1,nbr_t - floor(nbr_t/2))]);

% movie struct
clearvars F;
F(nbr_t) = struct('cdata',[],'colormap',[]);

figure(1); clf;
pos = get(gcf,'Position');
u2 = reshape(u,J1,J2);
u2(u2 < 0.01) = nan;
surf(X,Y,u2,'EdgeColor','none');
hold on
image([y(1), y(end)],[x(1), x(end)], country);
axis([y(1), y(end), x(1), x(end), 0, 10])
shading interp
axis ij
camorbit(dtheta(1),dphi(1))
camzoom(zo(1))
camdolly(0,dy(1),0)
axis off
drawnow;
F(1) = getframe(gcf,[0,0,pos(3:4)]);
disp('press any key to continue');
pause

% Crank-Nicolson implicit scheme 
A = (speye(J) - dt/h^2*D/2*L);

% Neumann no flux conditions: compute (approximate) boundary normal vector 
normal = zeros(size(border));
for i = 1:length(border)
    if any(interior == border(i)+1) % normal up ^
        normal(i) = 1;
    elseif any(interior == border(i)-1) % normal down v
        normal(i) = 3;
    elseif any(interior == border(i)-J1) % normal right >
        normal(i) = 2;
    elseif any(interior == border(i)+J1) % normal left <
        normal(i) = 4;
    else
        normal(i) = -1;
    end
end

vi = 2;

% BOUCLE PRINCIPALE
tic
while t < tfinal
    newu =A\(u + dt/h^2*D/2*L*u);
    
    % Neumann no flux
    for i = 1:length(border)
       if ( normal(i) == 1 )
           newu(border(i)) = newu(border(i)+1);
       elseif ( normal(i) == 2 )
           newu(border(i)) = newu(border(i)-J1);
       elseif ( normal(i) == 3 )
           newu(border(i)) = newu(border(i)-1);
       elseif ( normal(i) == 4 )
           newu(border(i)) = newu(border(i)+J1);
       else
           newu(border(i)) = nan;
       end     
    end
    
    u = newu;

    hold off
    u2 = reshape(u,J1,J2);
    u2(u2 < 0.01) = nan;
    surf(X,Y,u2,'EdgeColor','none');
    hold on
    image([y(1), y(end)],[x(1), x(end)], country);
    axis([y(1), y(end), x(1), x(end), 0, 10])
    shading interp
    axis ij
    axis off
    camorbit(dtheta(vi),dphi(vi))
    camzoom(zo(vi));
    camdolly(0,dy(vi),0);
    drawnow;
    F(vi) = getframe(gcf,[0,0,pos(3:4)]);
    t = t + dt;
    vi = vi + 1;
    fprintf("t = %.5f\n",t);
end
toc

%% Export simulation as video

v = VideoWriter('epidemics','MPEG-4');
v.Quality = 25;
open(v);
writeVideo(v,F);
close(v);

