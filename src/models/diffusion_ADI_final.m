% Heat equations finite difference using ADI
% Domain of integration: Metropolitan France (datasets/grid.mat)
% Miléna Kaag & Vincent Le Goff
%
% du/dt = D (d^2 u/dx^2 + d^2 u/dy^2)
%
% With Neumann no flux boundary conditions
% Initial conditions: 4 local peaks

% Equation parameters 
D = 400; % diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))

% Spatial grid in grid.mat

load ../datasets/grid.mat
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
dt = 0.2;

% camera zoom and orbit
nbr_t = ceil((tfinal-t0)/dt)+2; 
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

% ADI's Crank-Nicolson implicit schemes 
Ax=(speye(J)- dt/2/h^2*D*Lx);
Ay=(speye(J)- dt/2/h^2*D*Ly);

% SB: mldivide '\' ne semble par marcher, alors on va faire 
% une pre-factorisation LU: Ax = llx * uux et Ay = lly * uuy
% avec ll et uu des matrices triangulaires (similaire a la factorisation
% QR, mais legeremnet plus rapide
% le systeme lineaire A*u = b se resoud u = L\(U\b)
[llx,uux] = lu(Ax); 
[lly,uuy] = lu(Ay);



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

% MAIN LOOP
tic
while t < tfinal
    
    % SB: on a remplace les resolutions du systeme lineaire Ax\b par 
    % une factorisation manuelle LU: uux\(llx\b)
    % et Ay\b par uuy\(lly\b)
    b = ( u + dt/2/h^2*D*Ly*u ); 
    
    %Convert data in X dominant format
    b=reshape(reshape(b,J1,J2)',J,1); 
    u12=uux\(llx\b); 
    Lxu12 = Lx*u12; 
    
    %Convert data in Y dominant format
    u12 = reshape(reshape(u12,J2,J1)',J,1);
    Lxu12 = reshape(reshape(Lxu12,J2,J1)',J,1);
    newu = uuy\(lly\( u12 + dt/2/h^2*D*Lxu12));
    
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
    t = t + dt;
    fprintf("t = %.5f vi=%d\n",t,vi);
    
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
    vi = vi + 1;
    
end
toc

%% Export simulation as video

v = VideoWriter('epidemics','MPEG-4');
v.Quality = 25;
open(v);
writeVideo(v,F);
close(v);

