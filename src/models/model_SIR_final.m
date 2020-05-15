function sol = model_SIR()
% MODEL_SIR simulation of the SIR model with n classes
% MaeToday at 9:13 AM
% Pour le modele SIR  on a ça@Pauline , @Alexandre Vila , @Nino , @Ian


% Dynamical parameters
N = [2.0e8 ; 3.0e8 ; 1.0e8];    % population
%mu = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
%Lambda = N*mu;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
gam = [1/14.0 ; 1/15.0 ; 1/10.0];  % recovery rate 1/two weeks = 1/14 
beta = [0.17 0.3 0.1 ;
        0.50 0.2 0.21 ;
        0.40 0.12 0.23 ];   % infection rate S -> I

% Integration parameters 
IC = [(N(1)-0.25e4) 0.25e4 0 ; (N(2)-0.33e4) 0.33e4 0 ; (N(3)-0.42e4) 0.42e4 0 ]; % seed a few infected individuals
tspan = [0 365]; % in days 

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options);
soly = [sum(sol.y(1:3,:)) ;
        sum(sol.y(4:6,:)) ;
        sum(sol.y(7:9,:)) ];

figure(2); clf;
semilogy(sol.x, soly(1:3,:));
legend('S','I','R');

disp(['Percentage of contaminated 100*(E+I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(2:3,end)/N))]);


    function dxdt = sir(~,xx)  % nested function 
        % SEIR is the ODE rhs 
        S = xx(1:3);
        I = xx(4:6); 
        R = xx(7:9);
        
        
        dxdt = [ (- beta * I .* S./N);
                 (beta * I .* S./N - gam .* I);
                 (gam .* I) ]; 
             
    end                         % end nested function sir

end                             % end main function run_SIR 
