% Heat equations finite difference implicit Crank-Nicolson

% parametre des equations
D = 10; % coefficient de diffusion

% discretisation du Laplacien avec conditions aux bord de Neumann
% no flux, i.e. du/dx = 0 en x0 et en x1

load ../datasets/grid.mat
load ../laplacian/L2D.mat

% parametres de simulation, espace
% x in [x0,x1]
h = 1; % not relevant at the moment
x = 0:h:J2-1; % discretisation
y = 0:h:J1-1; 
[X,Y] = meshgrid(x,y);
J  = J1*J2;  % nombre de points de discretisation

% variable dynamiques
u = zeros(J,1); % stocke seulement l'etat au temps t
newu = zeros(J,1);



% condition initiales
u( (X(:) - 50).^2 + (Y(:) - 50).^2 < 3 ) = 1; % condition initiale constante par morceaux

exterior = setdiff(1:J, union(border,interior));

u(exterior) = nan;

X(exterior) = nan;
Y(exterior) = nan;

% parametres de simulation, temps
t0 = 0;
tfinal = 30; 
t = t0;
tp = t0;
dt = 0.5 * h^2/2/D;

nbr_t = ceil((tfinal-t0)/dt)+1;
dtheta = linspace(-37.5,37.5,floor(nbr_t/2));
dphi = linspace(-10,60,floor(nbr_t/2));
zo = linspace(5,1,floor(nbr_t/2));
dy = linspace(-3,0,floor(nbr_t/2));
dtheta = [dtheta, repmat(37.5,1,nbr_t - floor(nbr_t/2))];
dphi = [dphi, repmat(60,1,nbr_t - floor(nbr_t/2))];
zo = [zo, repmat(1,1,nbr_t - floor(nbr_t/2))];
dy = [dy, repmat(0,1,nbr_t - floor(nbr_t/2))];

clearvars F;
F(nbr_t) = struct('cdata',[],'colormap',[]);

figure(1); clf;
surf(X,Y,reshape(u,J1,J2),'EdgeColor','none');
axis([1, J1, 1, J2, 0, 10])
shading interp
axis ij
camorbit(dtheta(1),dphi(1))
camzoom(zo(1))
camdolly(0,dy(1),0)
axis off
drawnow;
F(1) = getframe(gcf,[0,0,853,544]);
disp('appuyer sur une touche pour continuer');
pause

% Schema implicite Crank-Nicolson
A = (speye(J) - dt/h^2*D/2*L);

% Neumann no flux conditions
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

    surf(X,Y,reshape(u,J1,J2),'EdgeColor','none');
    axis([1, J1, 1, J2, 0, 10])
    shading interp
    axis ij
    axis off
    camorbit(dtheta(vi),dphi(vi))
    camzoom(zo(vi));
    camdolly(0,dy(vi),0);
    drawnow;
    F(vi) = getframe(gcf,[0,0,853,544]);
    t = t + dt;
    vi = vi + 1;
    fprintf("t = %.5f\n",t);
end
toc

%%

v = VideoWriter('epidemics.avi','MPEG-4');
v.Quality = 25;
open(v);
writeVideo(v,F);
close(v);

