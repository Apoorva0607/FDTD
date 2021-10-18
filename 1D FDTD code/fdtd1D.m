% this code is to test 1D FDTD simulation with a square pulse on from 1 to
% 20 time steps
clc;
clear;
%% Input S (courant stability factor)
n = input('Enter the parameter S(courant stabilityb factor): ');

switch n
    case 1
        S = 0.5;
    case 2
        S = 1.01;
    case 3
        S = 0.99;
    
end

%% Initializing variables
eps = 8.85418782*1e-12; % m-3 kg-1 s4 A2 permittivity of free space
mu = 1.25663706*1e-6; % m kg s-2 A-2 permeability of free space
Da = 1; % constant from faraday's law
Ca = 1; % constant from ampere's law
delta = 1; % dx
c = 2.99792458*1e+08;% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
Db = dt/(mu*delta); % constant from faraday's law
Cb = dt/(eps*delta); % constant from ampere's law
imax = 200;
nmax = 50;
Ez = zeros(imax,1); % Initializing the array
Hy = zeros(imax,1); % Initializing the array
%%
for N = 1:nmax
    % defining the source function
    if N <= 20
        Ez(imax/4) = 1;
    else
        Ez(imax/4) = 0;
    end
    
    
    for i = 1:imax-1  % update H field
        Hy(i) = Da*Hy(i) + Db*(Ez(i+1)-Ez(i));
        
    end
    Ez(imax) = 0;
    for i = 1:imax-1  % update E field
        
        Ez(i+1) = Ca*Ez(i+1) + Cb*(Hy(i+1)-Hy(i));
    end
    
    
    % plotting Ez 
    figure(2);
    plot(Ez,'LineWidth',2);
    set(gca, 'FontSize',14);
    title('Free space grid', 'FontSize', 14);
    ylabel('Ez(i)', 'FontSize', 14);
    xlabel('Grid i coordinate','FontSize', 14);
    
end

