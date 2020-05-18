%% COVID-19 epidemiological reaction-diffusion-jump model
% Reaction: structured SIR model
% Diffusion: Metropolitan France, 100km^2 resolution
% Jump: airports/train stations/highway trafic
% 
% Numerical schemes
% Reaction: Forward Euler 
% Diffusion: Implicit Crank-Nicolson finite difference using 
% alternating direction implicit method (ADI)
% Domain of integration: Metropolitan France (datasets/grid.mat)
%                        Neumann no flux boundary conditions 
% Transport: jumps from and to transport hubs (train stations, airport and
% city centers) 
%
% Contribution:
% ADI scheme: Milena Kaag & Vincent Le Goff
% SIR model: @Pauline , @Alexandre Vila , @Nino , @Ian
% Diffusion: @Michel @Raf @Julie @Adrian et @Loup
% Transport: Mathieu, Virgile, flo, Mallou, Alexandre Ve
% du/dt = D (d^2 u/dx^2 + d^2 u/dy^2)

%% Diffusion parameters 
% SB: Ici on va faire des hypothèses pratiques mais pas vraiment réalistes:
%     - Les susceptibles et les removed ne diffusent pas du tout
%     - Les infectés diffusent (et peuvent donc aller infecter des populations éloignées)
%     La raison pour faire ces hypothèses est qu'on va supposer une densité inhomogène de 
%    susceptibles sur le territoire (Ville vs regions rurales) et qu'on veut éviter une exode des 
%    villes vers la campagne juste a cause de la diffusion. Tant que le nombre d'infectes est petit
%    la densité totale ne sera pas trop affectée.
%    Pour plus de realisme on pourra utiliser les reseaux de transport mis en place par le groupe transport
DS = 0;
DI = 8;
DR = 0;

%% SIR model parameters
gam = [1/14.0 ; 1/15.0 ; 1/10.0];  % recovery rate 1/two weeks = 1/14 
Gam = diag(gam);
beta = [0.17 0.3 0.1 ;
        0.50 0.2 0.21 ;
        0.40 0.12 0.23 ];   % infection rate S -> I


%% Spatial discretization
load ../datasets/grid.mat
load ../datasets/country.mat
load ../laplacian/L2D_ADI.mat

%% Simulation parameters 
% data in datasets/grid.mat
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
ymax = 3000;

exterior = setdiff(1:J, union(border,interior));

X(exterior) = nan;
Y(exterior) = nan;


%% Dynamical variables 


%% Initial Conditions
basal_populations_density = 60; % Default population density in France 
basal_infection_density = 0.0;

% Grandes villes
% columns: X pos; Y pos; Radius; population density, infected densities
cities = [
550, 250, 30.09, 3676.9, 0.1;  % Paris @ 3676.9 hab/km^2
753, 615, 30.09, 1317, 0.001;    % Lyon @ 1317 hab/km^2 
815, 850, 23.47, 900.6,  0.0;   % Marseille @ 900.6 hab/km^2 
472, 855, 16.07, 1083.9, 0.0;  % Toulouse @ 1083.9 hab/km^2
320, 659, 19.32, 719.2, 0.01;   % Bordeaux @ 719.2 hab/km^2
600, 50 , 11.87, 2301.3, 0.0;  % Lille @ 2301.3 hab/km^2
925, 815, 15.38, 1266.5, 0.0;  % Nice @ 1266.5 hab/km^2
220, 425, 13.08, 1100,   0.0;    % Nantes @ 1100 hab/km^2
945, 266, 8.74 , 1873.5, 2.0;  % Strasbourg @ 1873.5 hab/km^2
240, 330, 9.52 , 1092.4, 0.1;  % Rennes @ 1092.4 hab/km^2
831, 663, 12.77, 970,    0.0];    % Grenoble @ 970 hab/km^2

%Vide
empty = [
652,  433, 80, 1.0, 0.0;  % Entre Lyon et Paris @ 1 hab/km^2
580,  700, 80, 1.0, 0.0]; % Massif Central @ 1 hab/km^2

% dynamical variables, initial values
sys_size = 9;
nclass = [2.0e8 ; 3.0e8 ; 1.0e8];    % population distribution
nclass = nclass/sum(nclass); 
% S(1:3), I(1:3), R(1:3)
u = zeros(J,sys_size);
newu = zeros(J,sys_size);
dudt = zeros(J,sys_size);

% Fill regions and cities with susceptible and infected
u(:,1:3) = basal_populations_density;
u(:,4:6) = basal_infection_density;
for i=1:size(cities,1)
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 1 ) = cities(i,4) * nclass(1); 
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 2 ) = cities(i,4) * nclass(2); 
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 3 ) = cities(i,4) * nclass(3); 
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 4 ) = cities(i,5) * nclass(1); 
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 5 ) = cities(i,5) * nclass(2); 
  u( (X(:) - cities(i,1)).^2 + (Y(:) - cities(i,2)).^2 < cities(i,3)^2, 6 ) = cities(i,5) * nclass(3); 
end
for i=1:size(empty,1)
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 1 ) = empty(i,4) * nclass(1); 
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 2 ) = empty(i,4) * nclass(2); 
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 3 ) = empty(i,4) * nclass(3); 
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 4 ) = empty(i,5) * nclass(1); 
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 5 ) = empty(i,5) * nclass(2); 
  u( (X(:) - empty(i,1)).^2 + (Y(:) - empty(i,2)).^2 < empty(i,3)^2, 6 ) = empty(i,5) * nclass(3); 
end

u(exterior,:) = nan;

% Diffusion coefficient matrix
D = diag([DS; DS; DS; DI; DI; DI; DR; DR; DR]);

%% time parameters
t0     = 0;
tfinal = 90; 
t      = t0;
dt     = 1;

%% Display 
% camera zoom and orbit
nbr_t = ceil((tfinal-t0)/dt)+2; 
% dtheta = linspace(-37.5,37.5,floor(nbr_t/2));
% dphi = linspace(-10,60,floor(nbr_t/2));
% zo = linspace(5,1,floor(nbr_t/2));
% dy = linspace(-3,0,floor(nbr_t/2));
% dtheta = fliplr([dtheta, repmat(37.5,1,nbr_t - floor(nbr_t/2))]);
% dphi = fliplr([dphi, repmat(60,1,nbr_t - floor(nbr_t/2))]);
% zo = fliplr([zo, ones(1,nbr_t - floor(nbr_t/2))]);
% dy = fliplr([dy, zeros(1,nbr_t - floor(nbr_t/2))]);

% movie struct
clearvars F;
F(nbr_t) = struct('cdata',[],'colormap',[]);

figure(1); clf;
pos = get(gcf,'Position');

% PLOT THE TOTAL DENSITY OF INFECTED: SUM(U(:,4:6))
subplot(3,6,[1,2,3,7,8,9,13,14,15]);
u_out = reshape(sum(u(:,4:6),2),J1,J2);
u_out(u_out < 0.01) = nan; % don't display density if < 0.01
u_out(1,1) = 1000;
image([y(1), y(end)],[x(1), x(end)], country);
hold on
contourf(X,Y,log10(u_out),[-2,-1,0,1,2,3],'LineW',0.5);
axis([y(1), y(end), x(1), x(end), 0, ymax])
shading interp
axis ij
% camorbit(dtheta(1),dphi(1))
% camzoom(zo(1))
% camdolly(0,dy(1),0)
axis off
view(2);
axis equal;
drawnow;
F(1) = getframe(gcf,[0,0,pos(3:4)]);

% total densities
tt = 0:dt:tfinal;
utot = zeros(length(tt),sys_size);
utot(1,:) = sum(u(interior,:),1)*h^2;
panels_utot = [4,5,6,10,11,12,16,17,18];
titles_out = {'S_1','S_2','S_3','I_1','I_2','I_3','R_1','R_2','R_3'};
for i=1:sys_size
    subplot(3,6,panels_utot(i));
    area(tt,utot(:,i));
    title(titles_out{i});
end


disp('press any key to continue');
pause

vi = 2;

%% ADI's Crank-Nicolson implicit schemes 
ASx = (speye(J) - dt/h^2*DS/2*Lx);
AIx = (speye(J) - dt/h^2*DI/2*Lx);
ARx = (speye(J) - dt/h^2*DR/2*Lx);
ASy = (speye(J) - dt/h^2*DS/2*Ly);
AIy = (speye(J) - dt/h^2*DI/2*Ly);
ARy = (speye(J) - dt/h^2*DR/2*Ly);

% SB: mldivide '\' ne semble par marcher, alors on va faire 
% une pre-factorisation LU: Ax = llx * uux et Ay = lly * uuy
% avec ll et uu des matrices triangulaires (similaire a la factorisation
% QR, mais legeremnet plus rapide
% le systeme lineaire A*u = b se resoud u = L\(U\b)
[llsx,uusx] = lu(ASx); 
[llix,uuix] = lu(AIx); 
[llrx,uurx] = lu(ARx); 
[llsy,uusy] = lu(ASy);
[lliy,uuiy] = lu(AIy);
[llry,uury] = lu(ARy);

%% Neumann no flux conditions: compute (approximate) boundary normal vector 
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


%% MAIN LOOP
tic
while t < tfinal
    
    % SB: on a remplace les resolutions du systeme lineaire Ax\b par 
    % une factorisation manuelle LU: uux\(llx\b)
    % et Ay\b par uuy\(lly\b)
      
    dudt(:,1:3) = - (u(:,4:6)*beta').*u(:,1:3)./(sum(u,2));
    dudt(:,4:6) =   (u(:,4:6)*beta').*u(:,1:3)./(sum(u,2)) - u(:,4:6) * Gam;
    dudt(:,7:9) =    u(:,4:6) * Gam;

    b = ( u + dt/2 * dudt + dt/2/h^2*Ly*u*D ); 
    
    %Convert data in X dominant format
    b=reshape(permute(reshape(b,J1,J2,sys_size),[2,1,3]),J,sys_size); 
    u12(:,1:3)=uusx\(llsx\b(:,1:3)); 
    u12(:,4:6)=uuix\(llix\b(:,4:6)); 
    u12(:,7:9)=uurx\(llrx\b(:,7:9)); 
    Lxu12 = Lx*u12; 
    
    %Convert data in Y dominant format
    u12 = reshape(permute(reshape(u12,J2,J1,sys_size),[2,1,3]),J,sys_size);
    Lxu12 = reshape(permute(reshape(Lxu12,J2,J1,sys_size),[2,1,3]),J,sys_size);
    dudt(:,1:3) = - (u12(:,4:6)*beta').*u12(:,1:3)./(sum(u12,2));
    dudt(:,4:6) =   (u12(:,4:6)*beta').*u12(:,1:3)./(sum(u12,2)) - u12(:,4:6) * Gam;
    dudt(:,7:9) =    u12(:,4:6) * Gam;
    b = u12 + dt/2 * dudt + dt/2/h^2*Lxu12*D;
    newu(:,1:3) = uusy\(llsy\b(:,1:3));
    newu(:,4:6) = uuiy\(lliy\b(:,4:6));
    newu(:,7:9) = uury\(llry\b(:,7:9));

    % Neumann no flux
    for i = 1:length(border)
       if ( normal(i) == 1 )
           newu(border(i),:) = newu(border(i)+1,:);
       elseif ( normal(i) == 2 )
           newu(border(i),:) = newu(border(i)-J1,:);
       elseif ( normal(i) == 3 )
           newu(border(i),:) = newu(border(i)-1,:);
       elseif ( normal(i) == 4 )
           newu(border(i),:) = newu(border(i)+J1,:);
       else
           newu(border(i),:) = nan;
       end     
    end
    
    u = newu;
    t = t + dt;
    fprintf("t = %.5f, Infected=%f\n",t,sum(sum(u(interior,4:6),2))*h^2);
    
    subplot(3,6,[1,2,3,7,8,9,13,14,15]);
    hold off
    u_out = reshape(sum(u(:,4:6),2),J1,J2);
    u_out(u_out < 0.01) = nan;
    u_out(1,1) = 1000;
    image([y(1), y(end)],[x(1), x(end)], country);
    hold on
    axis([y(1), y(end), x(1), x(end), 0, ymax])
    contourf(X,Y,log10(u_out),[-2,-1,0,1,2,3],'LineW',0.5);
    colorbar
    shading interp
    axis ij
    axis off
%     camorbit(dtheta(vi),dphi(vi))
%     camzoom(zo(vi));
%     camdolly(0,dy(vi),0);
    view(2);
    axis equal;
    title(['day ' num2str(t)]);
    drawnow;
    F(vi) = getframe(gcf,[0,0,pos(3:4)]);
    
    % total densities
    utot(vi,:) = sum(u(interior,:),1)*h^2;
    for i=1:sys_size
        subplot(3,6,panels_utot(i));
        area(tt,utot(:,i));
        title(titles_out{i});
    end
    vi = vi + 1;
    
end
toc

%% Export simulation as video

F = F(1:vi-1);
v = VideoWriter('epidemics','MPEG-4');
v.Quality = 25;
open(v);
writeVideo(v,F);
close(v);

