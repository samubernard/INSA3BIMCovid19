%modèle SIR classique

% Paramètres du modèle 
beta = [0 2.1*10^-8 3*10^-8 ]; % taux d'infection
gamma = [0 1/7 1/12]; % taux de rémission
N = 85.7*10^5; % Population totale N = S + I + R
I0 = 791; % Nombre initial de personnes infectées
T1=0; %période préconfinement
T2=90; % période pendant le confinement
T3=200; %période post confinement
T=[T1 T2 T3] ; %temps total d'étude
dt = 0.0001; % pas de temps 

% appel à la fonction SIR_model
[S,I,R] = sir_model(beta,gamma,N,I0,T,dt); 

% Affichage du modèle

i = 0:dt:sum(T)-dt; % axe des X avec un pas de temps de dt et une valeur finale T-dt
plot(i,S,'b',i,I,'r',i,R,'g','LineWidth',2);
%ylim([0;10000])
grid on;
xlabel('Jours'); 
ylabel('Nombre individus');
legend('S','I','R');


%tableau avec les données importantes
% Data=zeros(3,3);
% Data(1,1)=I(1);
% Data(1,2)=R(1);
% Data(1,3)=I(1)+R(1);
% Data(2,1)=I(T1/dt);
% Data(2,2)=R(T1/dt);
% Data(2,3)=I(T1/dt)+R(T1/dt);
% Data(3,1)=I(T/dt);
% Data(3,2)=R(T/dt);
% Data(3,3)=I(T/dt)+R(T/dt);


% Calcul des constantes associées au modèle (R0 etc...)

% fonction SIR

function [S,I,R] = sir_model(beta,gamma,N,I0,T,dt)
    
    % Création d'un vecteur population susceptible
    S = zeros(1,sum(T)/dt); % crée un vecteur de 0 comportant 1 ligne et T/dt colonnes
    S(1) = N; % Le 1er élément du vecteur vaut la population : initialement, toute la population est susceptible
    % Création d'un vecteur population infectée
    I = zeros(1,sum(T)/dt); % crée un vecteur de 0 comportant 1 ligne et T/dt colonnes
    I(1) = I0; % Le 1er élément du vecteur vaut le nombre initial d'infectés : initialement, seulement I0 personnes sont infectées
    % Création d'un vecteur population guéries
    R = zeros(1,sum(T)/dt);% crée un vecteur de 0 comportant 1 ligne et T/dt colonnes
    % au début de l'étude, il n'y a pas de personnes guéries donc R(1)=0
    for i= 1:(sum(T)/dt)-1 %T/dt-1 itérations
        % Equations du modèle
        if i<T(1,1)/dt; 
            beta_util=beta(1,1);
            gamma_util=gamma(1,1);
        else if i>=T(1,1)/dt && i<((T(1,1)+T(1,2))/dt);
                beta_util=beta(1,2);
                gamma_util=gamma(1,2);
            else i>=(T(1,1)+T(1,2))/dt && i<sum(T)/dt;
                beta_util=beta(1,3);
                gamma_util=gamma(1,3);
            end
        end
        dS = (-beta_util*I(i)*S(i)) * dt;
        dI = (beta_util*I(i)*S(i) - gamma_util*I(i)) * dt;
        dR = (gamma_util*I(i)) * dt;
        S(i+1) = S(i) + dS;
        I(i+1) = I(i) + dI;
        R(i+1) = R(i) + dR;
    end

end
