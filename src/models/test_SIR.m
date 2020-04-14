% Heat equations finite difference implicit Crank-Nicolson
% Domain of integration: Metropolitan France (datasets/grid.mat)
% Discretisation: laplacian/L2D.math
%
% du/dt = D (d^2 u/dx^2 + d^2 u/dy^2)
%
% With Neumann no flux boundary conditions
% Initial conditions: 4 local peaks

% Equation parameters 
DS = 400; % Susceptible diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))
DI = 100; % Infected diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))
DR = 400; % Removed

% Dynamical parameters 
N = 6.0e8;      % total population
mu = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
Lambda = N*mu;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
a = 1/5.2;     % incubation rate E -> I 1/incubatio period = 1/5.2
gam = 1/14.0;  % recovery rate 1/two weeks = 1/14 

% transmission rate based on R0 = 2.5
% R0 = a*beta/(mu+a)/(mu+gam)
% beta = R0*(mu+a)*(mu+gam)/a
R0 = 2.5;
beta = R0*(mu+a)*(mu+gam)/a;   % infection rate S -> E


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
% S, I, R
S = zeros(J,1); 
I = zeros(J,1); 
R = zeros(J,1);
newS = zeros(J,1); 
newI = zeros(J,1);
newR = zeros(J,1);

% Initial Conditions
% Mostly Susceptible, and a few infected

I( (X(:) - 550).^2 + (Y(:) - 250).^2 < 300 ) = 1; % piecewise constant initial condition
I( (X(:) - 850).^2 + (Y(:) - 270).^2 < 300 ) = 1; % piecewise constant initial condition
I( (X(:) - 670).^2 + (Y(:) - 650).^2 < 300 ) = 1; % piecewise constant initial condition
I( (X(:) - 330).^2 + (Y(:) - 700).^2 < 300 ) = 1; % piecewise constant initial condition
S = N - I;
% R = 0 already;


exterior = setdiff(1:J, union(border,interior));

S(exterior) = nan;
I(exterior) = nan;
R(exterior) = nan;

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

% Plot infected
figure(1); clf;
pos = get(gcf,'Position');
I2 = reshape(I,J1,J2);
I2(I2 < 0.01) = nan;
surf(X,Y,I2,'EdgeColor','none');
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
AS = (speye(J) - dt/h^2*DS/2*L);
AI = (speye(J) - dt/h^2*DI/2*L);
AR = (speye(J) - dt/h^2*DR/2*L);

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

% SEIR model
%         Lambda - mu*S - beta*I/N*S;
%         beta*I/N*S - (mu+gam)*I;
%         gam*I - mu*R ]; 


% BOUCLE PRINCIPALE
tic
while t < tfinal
    ReactionS = Lambda - mu*S - beta*I.*S/N;
    ReactionI = beta*I.*S/N - (mu+gam)*I;
    ReactionR = gam*I - mu*R;
    newS =AS\(S + dt*ReactionS + dt/h^2*DS/2*L*S);
    newI =AI\(I + dt*ReactionI + dt/h^2*DI/2*L*I);
    newR =AR\(R + dt*ReactionR + dt/h^2*DR/2*L*R);
    
    % Neumann no flux
    for i = 1:length(border)
       if ( normal(i) == 1 )
           newS(border(i)) = newS(border(i)+1);
           newI(border(i)) = newI(border(i)+1);
           newR(border(i)) = newR(border(i)+1);
       elseif ( normal(i) == 2 )
           newS(border(i)) = newS(border(i)-J1);
           newI(border(i)) = newI(border(i)-J1);
           newR(border(i)) = newR(border(i)-J1);
       elseif ( normal(i) == 3 )
           newS(border(i)) = newS(border(i)-1);
           newI(border(i)) = newI(border(i)-1);
           newR(border(i)) = newR(border(i)-1);
       elseif ( normal(i) == 4 )
           newS(border(i)) = newS(border(i)+J1);
           newI(border(i)) = newI(border(i)+J1);
           newR(border(i)) = newR(border(i)+J1);
       else
           newS(border(i)) = nan;
           newI(border(i)) = nan;
           newR(border(i)) = nan;
       end     
    end
    
    S = newS;
    I = newI;
    R = newR;
    

    hold off
    I2 = reshape(I,J1,J2);
    I2(I2 < 0.01) = nan;
    surf(X,Y,I2,'EdgeColor','none');
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

