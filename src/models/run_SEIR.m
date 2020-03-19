function sol = run_SEIR()
% RUN_SEIR simulation of the SEIR model

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

% Integration parameters 
IC = [(N-1e4); 1e4; 0; 0]; % seed a few exposed individuals
tspan = [0 365]; % in days 

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@seir,tspan,IC,options);

figure(2); clf;
semilogy(sol.x, sol.y(1:4,:));

disp(['Percentage of contaminated 100*(E+I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(2:4,end)/N))]);


    function dxdt = seir(~,xx)  % nested function 
        % SEIR is the ODE rhs 
        
        S = xx(1); E = xx(2); I = xx(3); R = xx(4);
        
        dxdt = [ Lambda - mu*S - beta*I/N*S;
                 beta*I/N*S - (mu+a)*E;
                 a*E - (gam+mu)*I;
                 gam*I - mu*R ]; 
    end                         % end nested function seir

end                             % end main function run_SEIR 
