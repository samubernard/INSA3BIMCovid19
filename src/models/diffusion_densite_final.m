% GROUPE DENSITE
% MichelToday at 9:06 AM
% Pour la densit�, je suis avec @Raf @Julie @Adrian
% et @Loup

% SB: VOIR COMMENTAIRES CI-BAS

% Heat equations finite difference implicit Crank-Nicolson
% Domain of integration: Metropolitan France (datasets/grid.mat)
% Discretisation: laplacian/L2D.math
%
% du/dt = D (d^2 u/dx^2 + d^2 u/dy^2)
%
% With Neumann no flux boundary conditions
% Initial conditions: 4 local peaks

% Equation parameters 
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
% DS = 400; % Susceptible diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))
% DI = 100; % Infected diffusion coefficient units: km^2/day (characteristic travel speed sigma = sqrt(2*D))
% DR = 400; % Removed

% Dynamical parameters 
% SB: les unités de S,I et R sont en densité de population: habitant/km^2
%     il faut donc repartir les populations sur le territoire d'intérêt
%     Pour Paris (Ile de France), la superficie du territoire est 12000 k^2
%     avec 12 000 000 d'habitant => 1000 habitants/km^2
%     Pour le Grand Lyon, on a 533 km^2 avec 1 400 000 habitants => 2626 hab/km^2
%     si on suppose que le reste de la france a une
%     densité hors agglomeration de disons 60 habitants/km^2
%     on obtient une population initiale de susceptibles inegalement reperatie sur le territoire



% N1 = 30000000;
% N2 = 30000000;
% N = N1 + N2;      % total population



mu = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
% SB: Lambda va dépendre de l'espace, voir plus bas
% Lambda = N*mu;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
a = 1/5.2;     % incubation rate E -> I 1/incubatio period = 1/5.2
gam = 1/14.0;  % recovery rate 1/two weeks = 1/14 

% SB: pas nécessaire, une seule population suffit
% SIR2 : city 2
% Parameters for Lyon
% mu2 = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
% Lambda2 = N*mu2;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
% a2 = 1/5.2;     % incubation rate E -> I 1/incubatio period = 1/5.2
% gam2 = 1/14.0;  % recovery rate 1/two weeks = 1/14 

% transmission rate based on R0 = 2.5
% R0 = a*beta/(mu+a)/(mu+gam)
% beta = R0*(mu+a)*(mu+gam)/a
R0 = 2.5;
beta = R0*(mu+a)*(mu+gam)/a;   % infection rate S -> E

% parameters for city 2
% R02 = 2.5;
% beta2 = R02*(mu2+a2)*(mu2+gam2)/a2;


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
% SB: Ile de France 12 000 km^2 soit un 'disque' de 61.8038 de rayon 
%     Lyon: 533 km^2 soit un disque de 13.025313 de rayon
%     Note: pas besoin de définir deux populations différentes pour Paris et Lyon,
%     un seul système SIR est suffisant


% S, I, R pour Paris
% SB: S, I, R pour Paris ET Lyon
S = zeros(J,1); 
I = zeros(J,1); 
R = zeros(J,1);
newS = zeros(J,1); 
newI = zeros(J,1);
newR = zeros(J,1);

% SIR2
% SB: Pas nécessaire
% S2 = zeros(J,1); 
% I2 = zeros(J,1); 
% R2 = zeros(J,1);
% newS2 = zeros(J,1); 
% newI2 = zeros(J,1);
% newR2 = zeros(J,1);

% Initial Conditions
% Mostly Susceptible, and a few infected

% SB: Susceptibles
S(:) = 60; % densité par défaut sur la France

%Grandes villes
S( (X(:) - 550).^2 + (Y(:) - 250).^2 < 30.09^2 ) = 3676.9; % Paris @ 3676.9 hab/km^2
S( (X(:) - 753).^2 + (Y(:) - 615).^2 < 30.09^2 ) = 1317; % Lyon @ 1317 hab/km^2 
S( (X(:) - 815).^2 + (Y(:) - 850).^2 < 23.47^2 ) = 900.6; % Marseille @ 900.6 hab/km^2 
S( (X(:) - 472).^2 + (Y(:) - 855).^2 < 16.07^2 ) = 1083.9; % Toulouse @ 1083.9 hab/km^2
S( (X(:) - 320).^2 + (Y(:) - 659).^2 < 19.32^2 ) = 719.2; % Bordeaux @ 719.2 hab/km^2
S( (X(:) - 600).^2 + (Y(:) - 50).^2 < 11.87^2 ) = 2301.3; % Lille @ 2301.3 hab/km^2
S( (X(:) - 925).^2 + (Y(:) - 815).^2 < 15.38^2 ) = 1266.5; % Nice @ 1266.5 hab/km^2
S( (X(:) - 220).^2 + (Y(:) - 425).^2 < 13.08^2 ) = 1100; % Nantes @ 1100 hab/km^2
S( (X(:) - 945).^2 + (Y(:) - 266).^2 < 8.74^2 ) = 1873.5; % Strasbourg @ 1873.5 hab/km^2
S( (X(:) - 240).^2 + (Y(:) - 330).^2 < 9.52^2 ) = 1092.4; % Rennes @ 1092.4 hab/km^2
S( (X(:) - 831).^2 + (Y(:) - 663).^2 < 12.77^2 ) = 970; % Grenoble @ 970 hab/km^2

%Petits villages (mal plac�s)
%S( (X(:) - 613).^2 + (Y(:) - 343).^2 < 3.04^2 ) = 6.5; % Lucenay-le-duc @ 6.5 hab/km^2
%S( (X(:) - 795).^2 + (Y(:) - 665).^2 < 3.09^2 ) = 1.9; % Châteauneuf-d'Entraunes @ 1.9 hab/km^2
%S( (X(:) - 431).^2 + (Y(:) - 519).^2 < 4.62^2 ) = 5; % Tarnac @ 5 hab/km^2
%S( (X(:) - 282).^2 + (Y(:) - 653).^2 < 6.86^2 ) = 7.6; % Sore @ 7.6 hab/km^2
%S( (X(:) - 595).^2 + (Y(:) - 640).^2 < 2.99^2 ) = 6.5; % Pied de Borne @ 6.5 hab/km^2
%S( (X(:) - 643).^2 + (Y(:) - 142).^2 < 2.26^2 ) = 6.1; % Bayonville @ 6.1 hab/km^2
%S( (X(:) - 608).^2 + (Y(:) - 390).^2 < 4.07^2 ) = 14; % Anost @ 14 hab/km^2
%S( (X(:) - 88).^2 + (Y(:) - 246).^2 < 4.22^2 ) = 16; % Berrien @ 16 hab/km^2
%S( (X(:) - 453).^2 + (Y(:) - 802).^2 < 1.95^2 ) = 8.4; % Camurac @ 8.4 hab/km^2
%S( (X(:) - 739).^2 + (Y(:) - 587).^2 < 6.26^2 ) = 0.87; % St Christophe en Oisans @ 0.87 hab/km^2

%Vide
S( (X(:) - 652).^2 + (Y(:) - 433).^2 < 80^2 ) = 1; % Entre Lyon et Paris @ 1 hab/km^2
S( (X(:) - 580).^2 + (Y(:) - 700).^2 < 80^2 ) = 1; % Massif Central @ 1 hab/km^2

I(:) = 0; % pas d'infection hors Paris/Lyon 
I( (X(:) - 550).^2 + (Y(:) - 250).^2 < 30.09^2 ) = 0.1; % Paris central (100km^2) 
I( (X(:) - 753).^2 + (Y(:) - 615).^2 < 30.09^2 ) = 0.001; % Lyon central (100km^2) 
%I( (X(:) - 815).^2 + (Y(:) - 850).^2 < 23.47^2 ) = 0.1; % Marseille central (100km^2) 
%I( (X(:) - 472).^2 + (Y(:) - 855).^2 < 16.07^2 ) = 0.1; % Toulouse central (100km^2) 
I( (X(:) - 320).^2 + (Y(:) - 659).^2 < 19.32^2 ) = 0.01; % Bordeaux central (100km^2) 
%I( (X(:) - 600).^2 + (Y(:) - 50).^2 < 11.87^2 ) = 0.1; % Lille central (100km^2) 
%I( (X(:) - 925).^2 + (Y(:) - 815).^2 < 15.38^2 ) = 0.1; % Nice central (100km^2) 
%I( (X(:) - 220).^2 + (Y(:) - 425).^2 < 13.08^2 ) = 0.1; % Nantes central (100km^2) 
I( (X(:) - 945).^2 + (Y(:) - 266).^2 < 8.74^2 ) = 2.0; % Strasbourg central (100km^2) 
%I( (X(:) - 240).^2 + (Y(:) - 330).^2 < 9.52^2 ) = 0.1; % Rennes central (100km^2) 
%I( (X(:) - 831).^2 + (Y(:) - 663).^2 < 12.77^2 ) = 0.1; % Grenoble central (100km^2) 
% SB R = 0

% SB: Enfin on peut définir Lambda
Lambda = S*mu;  % birth rate: set to keep initial population S constant (ignoring deaths from covid-19)


% I( (X(:) - 550).^2 + (Y(:) - 250).^2 < 300 ) = 1; % paris piecewise constant initial condition
%I2( (X(:) - 600).^2 + (Y(:) - 80).^2 < 300 ) = 1; % lille piecewise constant initial condition
%I( (X(:) - 928).^2 + (Y(:) - 235).^2 < 300 ) = 10; % strasbourg piecewise constant initial condition
% I2( (X(:) - 753).^2 + (Y(:) - 615).^2 < 300 ) = 1; % lyon piecewise constant initial condition
% S = N1 - I;
% S2 = N2 - I2;
% R = 0 already;


exterior = setdiff(1:J, union(border,interior));

S(exterior) = nan;
I(exterior) = nan;
R(exterior) = nan;

% SIR 2 : Lille
S2(exterior) = nan;
I2(exterior) = nan;
R2(exterior) = nan;

X(exterior) = nan;
Y(exterior) = nan;

% time parameters
t0 = 0;
tfinal = 440; 
t = t0;
dt = 1.0;

% camera zoom and orbit
% nbr_t = ceil((t1ce(-37.5,37.5,floor(nbr_t/2));
% dphi = linspace(-10,60,floor(nbr_t/2));
% zo = linspace(5,1,floor(nbr_t/2));
% dy = linspace(-3,0,floor(nbr_t/2));
% dtheta = fliplr([dtheta, repmat(37.5,1,nbr_t - floor(nbr_t/2))]);
% dphi = fliplr([dphi, repmat(60,1,nbr_t - floor(nbr_t/2))]);
% zo = fliplr([zo, ones(1,nbr_t - floor(nbr_t/2))]);
% dy = fliplr([dy, zeros(1,nbr_t - floor(nbr_t/2))]);

% movie struct
% clearvars F;
% F(nbr_t) = struct('cdata',[],'colormap',[]);

% Plot infected
figure(1); clf;
pos = get(gcf,'Position');
% Ii = reshape(I+I2,J1,J2);
Ii = reshape(I,J1,J2);
Ii(Ii < 0.01) = nan; % SB: on ne montre pas la solution si Ii < 0.01 infecté par km^2
surf(X,Y,Ii,'EdgeColor','none');
hold on
image([y(1), y(end)],[x(1), x(end)], country);
axis([y(1), y(end), x(1), x(end), 0, 100])
shading interp
axis ij
% camorbit(dtheta(1),dphi(1))
% camzoom(zo(1))
% camdolly(0,dy(1),0)
axis off
drawnow;
% F(1) = getframe(gcf,[0,0,pos(3:4)]);
disp('press any key to continue');
pause

% Crank-Nicolson implicit scheme 
AS = (speye(J) - dt/h^2*DS/2*L);
AI = (speye(J) - dt/h^2*DI/2*L);
AR = (speye(J) - dt/h^2*DR/2*L);

% SB: Pas nécessaire
%SIR 2 % pour l'instant inchang� (car on ne change pas DS, DI, DR pour
%l'instant)
% AS2 = (speye(J) - dt/h^2*DS/2*L);
% AI2 = (speye(J) - dt/h^2*DI/2*L);
% AR2 = (speye(J) - dt/h^2*DR/2*L);

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

    % SB: La population va diffuser, la densité va donc évoluer et S+I+R ne sera pas constante
    ReactionS = Lambda - mu*S - beta*I.*S./(S+I+R); % ce qui �volue � chaque pas de temps
    ReactionI = beta*I.*S./(S+I+R) - (mu+gam)*I;
    ReactionR = gam*I - mu*R;
    newS =AS\(S + dt*ReactionS + dt/h^2*DS/2*L*S); % on prend en compte le temps et l'espace
    newI =AI\(I + dt*ReactionI + dt/h^2*DI/2*L*I);
    newR =AR\(R + dt*ReactionR + dt/h^2*DR/2*L*R);


    % ReactionS = Lambda - mu*S - beta*I.*S/N1; % ce qui �volue � chaque pas de temps
    % ReactionI = beta*I.*S/N1 - (mu+gam)*I;
    % ReactionR = gam*I - mu*R;
    % newS =AS\(S + dt*ReactionS + dt/h^2*DS/2*L*S); % on prend en compte le temps et l'espace
    % newI =AI\(I + dt*ReactionI + dt/h^2*DI/2*L*I);
    % newR =AR\(R + dt*ReactionR + dt/h^2*DR/2*L*R);
    
    % SB: pas nécessaire
    %SIR 2
    % ReactionS2 = Lambda2 - mu2*S2 - beta2*I2.*S2/N2; 
    % ReactionI2 = beta2*I2.*S2/N2 - (mu2+gam2)*I2;
    % ReactionR2 = gam2*I2 - mu2*R2;
    % newS2 =AS2\(S2 + dt*ReactionS2 + dt/h^2*DS/2*L*S2);
    % newI2 =AI2\(I2 + dt*ReactionI2 + dt/h^2*DI/2*L*I2);
    % newR2 =AR2\(R2 + dt*ReactionR2 + dt/h^2*DR/2*L*R2);
    
    % Neumann no flux
    for i = 1:length(border)
       if ( normal(i) == 1 )
           newS(border(i)) = newS(border(i)+1);
           newI(border(i)) = newI(border(i)+1);
           newR(border(i)) = newR(border(i)+1);
           
           % newS2(border(i)) = newS2(border(i)+1);
           % newI2(border(i)) = newI2(border(i)+1);
           % newR2(border(i)) = newR2(border(i)+1);
           
       elseif ( normal(i) == 2 )
           newS(border(i)) = newS(border(i)-J1);
           newI(border(i)) = newI(border(i)-J1);
           newR(border(i)) = newR(border(i)-J1);
           
           % newS2(border(i)) = newS2(border(i)-J1);
           % newI2(border(i)) = newI2(border(i)-J1);
           % newR2(border(i)) = newR2(border(i)-J1);
           
       elseif ( normal(i) == 3 )
           newS(border(i)) = newS(border(i)-1);
           newI(border(i)) = newI(border(i)-1);
           newR(border(i)) = newR(border(i)-1);
           
           % newS2(border(i)) = newS2(border(i)-1);
           % newI2(border(i)) = newI2(border(i)-1);
           % newR2(border(i)) = newR2(border(i)-1);
           
       elseif ( normal(i) == 4 )
           newS(border(i)) = newS(border(i)+J1);
           newI(border(i)) = newI(border(i)+J1);
           newR(border(i)) = newR(border(i)+J1);
           
           % newS2(border(i)) = newS2(border(i)+J1);
           % newI2(border(i)) = newI2(border(i)+J1);
           % newR2(border(i)) = newR2(border(i)+J1);
           
       else
           newS(border(i)) = nan;
           newI(border(i)) = nan;
           newR(border(i)) = nan;
           
           % newS2(border(i)) = nan;
           % newI2(border(i)) = nan;
           % newR2(border(i)) = nan;
           
       end     
    end
    
    S = newS;
    I = newI;
    R = newR;
    
    % S2 = newS2;
    % I2 = newI2;
    % R2 = newR2;
    
    % total S,I,R
    disp(h^2*[sum(S(interior)),sum(I(interior)),sum(R(interior))]);

    hold off
    % Ii = reshape(I+I2,J1,J2);
    Ii = reshape(I,J1,J2);
    Ii(Ii < 0.01) = nan;
    surf(X,Y,Ii,'EdgeColor','none');
    hold on
    image([y(1), y(end)],[x(1), x(end)], country);
    axis([y(1), y(end), x(1), x(end), 0, 100])
    shading interp
    axis ij
    axis off
    % camorbit(dtheta(vi),dphi(vi))
    % camzoom(zo(vi));
    % camdolly(0,dy(vi),0);
    drawnow;
    % F(vi) = getframe(gcf,[0,0,pos(3:4)]);
    t = t + dt;
    vi = vi + 1;
    fprintf("t = %.5f\n",t);
end
toc

%% Export simulation as video

% v = VideoWriter('epidemics','MPEG-4');
% v.Quality = 25;
% open(v);
% writeVideo(v,F);
% close(v);

