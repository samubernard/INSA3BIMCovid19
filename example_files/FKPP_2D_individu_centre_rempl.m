% FKPP - version individu-centre 2D avec remplacement
% Des individus sont distribues sur un intervalle
% Certains mutants sont porteurs d'un gene favorable qui donne
% un avantage de croissance r
% Tous les individus de deplacent selon une marche aleatoire avec 
% des sauts en sqrt(2*D) 
% Les wild-types (non-porteurs) on un taux de croissance 0
% A chaque reproduction d'un mutant, l'individu le plus
% proche est remplace par la progeniture
% Le nombre total d'individu est constant = N0
% la croissance est limitee par le fait la probabilite de remplacer un
% wild-type est de 1 - umut, ou umut est
% la densite de population de mutant
%
% prob(mutant en x se reproduise) = dt*r
% prob(remplacer un wild-type) ~ (1-u(x))
%
% NB: ksdensity pas utilise car trop lent!

%% parametre des equations
r = 0.5;        % taux de croissance
D = 0.1;        % coefficient de diffusion

%% parametres de simulation, espace
% x in [x0,x1]
h = 0.1;        % intervalle de discretisation spatiale
Sx = 20.0;      % longueur de l'intervalle en X
Sy = 15.0;      % longueur de l'intervalle en Y
x0 = 0;         % bord gauche de l'intervalle
x1 = Sx;         % bord droit de l'intervalle
xi = x0:h:x1;   % discretisation
JX = length(xi); % nombre de points de discretisation

y0 = 0;         % bord gauche de l'intervalle
y1 = Sy;         % bord droit de l'intervalle
yi = y0:h:y1;   % discretisation
JY = length(yi); % nombre de points de discretisation

%% Individu-centre
N0 = 30000;   % population totale

% chaque individu possede une position x et trait g
% g = 1 : gene favorable (mutants)
% g = 0 : pas de gene favorable
x = x0 + Sx*rand([N0,1]); % individus distribues sur [x0,x1]
y = y0 + Sy*rand([N0,1]); % individus distribues sur [y0,y1]
g = false(1,N0);
g( x < Sx/10 & y < Sy/10 ) = true; % individus sur la gauche de l'intervalle sont porteurs
w = zeros(1,N0);     % probabilite pour l'individu de se reproduire

%% parametres de simulation, temps
t0 = 0;
tfinal = 30.0; 
t = t0;
dt = 0.1

%% Initialisation
figure(1); clf;
% plot(x(g),y(g),'b.',x(~g),y(~g),'r.') % chaque point est un individu
% (porteur en bleu et wild-type en rouge)
plot(x(g),y(g),'b.') % chaque point est un individu porteur

%% BOUCLE PRINCIPALE
tic
while t < tfinal
    drawnow;
    % on calcule pour chaque mutant une probabilite de se reproduire
    % (les individus non porteurs ne se reproduisent pas)
    
        % on calcule pour chaque mutant une probabilite de se reproduire
    % (les individus non porteurs ne se reproduisent pas)
    w  = dt*r*g;                                   % probabilite de se reproduire
    irep = find(rand(1,N0) < w);                 % realisation 
    % un individu se reproduisant remplace l'individu le plus proche
    % (autre que lui-meme)
    % la probabilite de choisir un wild-type est environ 1-umut 
    for i = 1:length(irep) 
        [~,rempl] = min((x - x(irep(i))).^2 ...
           + (y - y(irep(i))).^2 ... 
           + 10*(x == x(irep(i))).*(y == y(irep(i))));
        g(rempl) = true;
    end
    
    % on deplace les individu en marche aleatoire
    x = x + sqrt(dt)*sqrt(2*D)*randn(size(x));
    x = abs(x); % condition en x0 reflechissant
    y = y + sqrt(dt)*sqrt(2*D)*randn(size(y));
    y = abs(y); % condition en x0 reflechissant
    
    % affichage
    % plot(x(g),y(g),'b.',x(~g),y(~g),'r.') % chaque point est un individu
    % (porteur en bleu et wild-type en rouge)
    plot(x(g),y(g),'b.') % chaque point est un individu porteur

    axis([0, Sx, 0, Sy]);
    t = t + dt;
    fprintf('t = %.5f, Nmut = %d\n',t, sum(g));
end
toc
