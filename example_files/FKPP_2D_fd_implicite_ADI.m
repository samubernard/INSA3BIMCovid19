% Equation FKPP difference finie 2D inplicite Crank-Nicolson avec
% approximation ADI: alternate direction implicit method

%% parametre des equations
r = 0.5; % taux de croissance
D = 0.1; % coefficient de diffusion

%% parametres de simulation, espace
% x in [x0,x1]
h = 0.025;  % intervalle de discretisation spatiale
Sx = 10.0;
Sy = 6.0;
x0 = 0;   % bord gauche de l'intervalle
x1 = Sx;  % bord droit de l'intervalle
y0 = 0;
y1 = Sy;
x = x0:h:x1; % discretisation
y = y0:h:y1; 
[X,Y] = meshgrid(x,y);
J1 = length(x); 
J2 = length(y);
J  = J1*J2;  % nombre de points de discretisation

%% variable dynamiques
u = zeros(J,1); % stocke seulement l'etat au temps t
u1 = zeros(J,1);

% discretisation du Laplacien avec conditions aux bord de Neumann
% no flux, i.e. du/dx = 0 en x0 et en x1

%% Le Laplacien discretise est une matrice L de taille JxJ
% En 1D elle est symetrique et tridiagonale

coinhautgauche = 1;
coinbasgauche = J2;
coinhautdroit = (J1-1)*J2+1;
coinbasdroit = J;
bordgauche = 2:J2-1;
bordhaut = J2+1:J2:J2*(J1-2)+1;
bordbas = 2*J2:J2:J2*(J1-1);
borddroit = J2*(J1-1)+2:J-1;
bord = [coinhautgauche, coinhautdroit, coinbasgauche, coinbasdroit, ...
    bordgauche, bordhaut, bordbas, borddroit];
interieurY = setdiff(1:J, bord);

[ix,iy] = ind2sub([J2,J1],interieurY);
interieurX = sub2ind([J1,J2],iy,ix);

%% test indices Y
% MY = zeros(J2,J1);
% MY(interieurY) = 1;
% spy(MY)
%% test indices X
% MX = zeros(J1,J2);
% MX(interieurX) = 1;
% spy(MX)

%% interieur
% On construit deux matrices du Laplacien discretise: LX et LY. LX est le Laplacien
% discretise en prenant les indices en X d'abord:
% 
%   1    2 3 4 5 ...   J1
%   J1+1 ...         2*J1
%   ...
%   J1*(J2-1)+1 ... J1*J2
%
% LY est le Laplacien discretise en prenant les indices le long de Y d'abord:
%  
%   1 J2+1 ...(J1-1)*J2+1
%   2 J2+2            ...
%   3 
%   ...               ...
%   J2 2*J2 3*J ... J1*J2
%
% LX et LY sont des matrices tridiagonales

LX = sparse(interieurX,interieurX,-2,J,J); % matrice creuse, compacte en memoire
LX = LX + sparse(interieurX,interieurX+1,1,J,J);
LX = LX + sparse(interieurX,interieurX-1,1,J,J);

LY = sparse(interieurY,interieurY,-2,J,J); % matrice creuse, compacte en memoire
LY = LY + sparse(interieurY,interieurY+1,1,J,J);
LY = LY + sparse(interieurY,interieurY-1,1,J,J);

%% condition initiales
% u( X(:).^2 + Y(:).^2 < 3 ) = 0.8; % condition initiale constante par
% morceaux
u = exp( - X(:).^2/0.5 - Y(:).^2/2 ); % condition initiale lisse


%% parametres de simulation, temps
t0 = 0;
tfinal = 20; 
t = t0;
tp = t0;
dt = .1;

%% Schema implicite Crank-Nicolson
% avec ADI: alternate direction implicit method
AX = (speye(J) - dt/h^2*D/2*LX);
AY = (speye(J) - dt/h^2*D/2*LY);

%% BOUCLE PRINCIPALE

figure(2); clf;
image(x,y,255*reshape(u,J2,J1));
disp('appuyer sur une touche pour continuer');
pause

tic
while t < tfinal
    drawnow;
    % Methode ADI:
    % On resout en 2 etapes: d'abord on fait par rapport a X:
    % u12 = u + dt/2*r*u*(1-u) + dt/2/h^2*D*(LX*u12 + LY*u)
    % ensuite par rapport a Y:
    % u = u12 + dt/2*r*u12*(1-u12) + dt/2/h^2*D*(LX*u12 + LY*u)
    % L'evaluation se fait a dt/2 a chaque etape
    
    b = ( u + dt/2*r*u.*(1-u) + dt/2/h^2*D*LY*u ); % En X: on calcule les donnee du systeme implicite AX*u12 = b
    b = reshape(reshape(b,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
    u12 = AX\b;                                    % on resout
    LXu12 = LX*u12;                                % On stocke LX*u12
    u12 = reshape(reshape(u12,J1,J2)',J,1);        % on reconvertit u12 en Y-dominant
    LXu12 = reshape(reshape(LXu12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
    u =   AY\( u12 + dt/2*r*u12.*(1-u12) + dt/2/h^2*D*LXu12 ); % En Y: on resout AY*u = b 
    u(bordgauche) = u(bordgauche+J2);              % condition Neumann no flux
    u(borddroit) = u(borddroit-J2);
    u(bordhaut) = u(bordhaut+1);
    u(bordbas) = u(bordbas-1);
    if ( fix(10*t) > fix(10*tp) )                  % affichage tout les 0.1 unite de temps 
        image(x,y,255*reshape(u,J2,J1));
    end
    tp = t;
    t = t + dt;
    fprintf("t = %.5f\n",t);
end
toc
