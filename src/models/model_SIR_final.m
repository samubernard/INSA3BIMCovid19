% modele SIR : Pauline , Maelys, Nino , Ian, Alexandre Vi

function sol = model_SIR()
% MODEL_SIR simulation of the SIR model with n classes

% Dynamical parameters
N = popAge([18 40 70]);    % population
%gam = [1/17.0 ; 1/23.0 ; 1/35.0 ; 1/30.0]; 
gam = ones(4,1)*1/10.0; % recovery rate 1/10days

%R0 = 2.3;       %Avant confinement
%R0 = 0.5;       %Pendant confinement
%R0 = 0.9;       %Apres confinement


Ci = [1 2 1.5 0 ; 2 4 3 1 ; 1.5 3 2.5 1 ; 0 1 1 0.5];        %contact possible avant confinement
Cc = [0. 1.5 1 0 ; 1.5 2 1 0 ; 1 1 1 0.5 ; 0 0 0.5 0.5];        %contact possible pendant confinement
Cpc = [0. 1.5 1 0 ; 1.5 3 1.5 0 ; 1 1.5 2 0.5 ; 0 0 0.5 0.5];        %contact possible apres confinement

tc = 20; %debut du confinement
tpc = 80; %fin du confinement

infect = 0.04; %infectiosite = Ro*gam/C 



% Integration parameters 
I0 = [0 ; 0 ; 1 ; 0];
IC = [(N(1)-I0(1)) I0(1) 0 ; (N(2)-I0(2)) I0(2) 0 ; (N(3)-I0(3)) I0(3) 0  ; (N(4)-I0(4)) I0(4) 0  ];    %infecte au debut

tspan = [0 365]; % in days 


% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options);
soly = [sum(sol.y(1:4,:)) ;
        sum(sol.y(5:8,:)) ;
        sum(sol.y(9:12,:));
        sum(sol.y)];

figure(1); clf;     %Modele
semilogy(sol.x, soly(1:4,:));
legend('S','I','R','Total');
title('Model');

solys = [sol.y(1,:) ;
        sol.y(2,:) ;
        sol.y(3,:);
        sol.y(4,:);
        sum(sol.y(1:4,:))];
    
figure(2); clf;     %Suceptible
semilogy(sol.x, solys(1:5,:));
legend('0-18ans','18-40ans','40-70ans','+70ans','Total');
title('Susceptible');

solyi = [sol.y(5,:) ;
        sol.y(6,:) ;
        sol.y(7,:);
        sol.y(8,:);
        sum(sol.y(5:8,:))];
    
figure(3); clf;     %Infectious
semilogy(sol.x, solyi(1:5,:));
legend('0-18ans','18-40ans','40-70ans','+70ans','Total');
title('Infectious');

solyr = [sol.y(9,:) ;
        sol.y(10,:) ;
        sol.y(11,:);
        sol.y(12,:);
        sum(sol.y(9:12,:))];
    
figure(4); clf;     %Retablis
semilogy(sol.x, solyr(1:5,:));
legend('0-18ans','18-40ans','40-70ans','+70ans','Total');
title('Retablis');


disp(['Percentage of contaminated 100*(I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(5:12,end))/sum(N))]);


    function dxdt = sir(t,xx)  % nested function 
        % SIR is the ODE rhs 
        S = xx(1:4);
        I = xx(5:8); 
        R = xx(9:12);
        beta = ((t<tc)*Ci+(t>=tc & t<tpc)*Cc+(t>=tpc)*Cpc)*infect;
        
        dxdt = [ (- beta * I .* S./N);
                 (beta * I .* S./N - gam .* I);
                 (gam .* I) ]; 
             
    end                         % end nested function sir

end                             % end main function run_SIR 