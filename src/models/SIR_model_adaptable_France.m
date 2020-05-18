%mod�le SIR classique

% Param�tres du mod�le 
beta = [3.1*10^-9 7*10^-10 0 ]; % taux d'infection
gamma = [1/14 1/14 0]; % taux de r�mission
N = 6.7*10^7; % Population totale N = S + I + R
I0 = 30000; % Nombre initial de personnes infect�es
T1=26; %p�riode pr�confinement
T2=65; % p�riode pendant le confinement
T3=0; %p�riode post confinement
T=[T1 T2 T3] ; %temps total d'�tude
dt = 0.0001; % pas de temps 

% appel � la fonction SIR_model
[S,I,R] = sir_model(beta,gamma,N,I0,T,dt); 

% Affichage du mod�le

i = 0:dt:sum(T)-dt; % axe des X avec un pas de temps de dt et une valeur finale T-dt

subplot(121);
%plot(i,S,'b',i,I,'r',i,R,'g','LineWidth',2);
plot(i,S,'b','LineWidth',2);
legend('Susceptibles');
subplot(122);
plot(i,I,'r',i,R,'g','LineWidth',2);
grid on;
xlabel('Jours'); 
ylabel('Nombre d`individus');
legend('Infect�s','R�tablis');




%tableau avec les donn�es importantes
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


% Calcul des constantes associ�es au mod�le (R0 etc...)

% fonction SIR

function [S,I,R] = sir_model(beta,gamma,N,I0,T,dt)
    
    % Cr�ation d'un vecteur population susceptible
    S = zeros(1,sum(T)/dt); % cr�e un vecteur de 0 comportant 1 ligne et T/dt colonnes
    S(1) = N; % Le 1er �l�ment du vecteur vaut la population : initialement, toute la population est susceptible
    % Cr�ation d'un vecteur population infect�e
    I = zeros(1,sum(T)/dt); % cr�e un vecteur de 0 comportant 1 ligne et T/dt colonnes
    I(1) = I0; % Le 1er �l�ment du vecteur vaut le nombre initial d'infect�s : initialement, seulement I0 personnes sont infect�es
    % Cr�ation d'un vecteur population gu�ries
    R = zeros(1,sum(T)/dt);% cr�e un vecteur de 0 comportant 1 ligne et T/dt colonnes
    % au d�but de l'�tude, il n'y a pas de personnes gu�ries donc R(1)=0
    for i= 1:(sum(T)/dt)-1 %T/dt-1 it�rations
        % Equations du mod�le
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
