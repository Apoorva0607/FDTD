% this code is to test 1D FDTD simulation with a square pulse and a
% guassian source pulse
clc;
clear;
clear all;
%% Input S (courant stability factor)
J = input('Enter the parameter S(courant stabilityb factor): ');
% ensures that numerical omega is real ( refer numerical dispersion
% relationship)

switch J
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
delta = 1*1e-3; % dx( for the avalanche case with ISM band 820Mhz - 980Mhz, with 10 grid/wavelength, wavelegnth = c/f_highest)
c = 2.99792458*1e+08;% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
Db = dt/(mu*delta); % constant from faraday's law
Cb = dt/(eps*delta); % constant from ampere's law
S1 = 1; % flag to choose square wave vs gaussian wave
if S1 ==1  % gaussian wave
    imax = 400; % need more grid point to compute field for gaussian wave
    nmax = 415; % found by plotting the gaussian function with time
    nhalf = 20;% ieally should be 2/(pi*26MHz) but this would make the spectrum go beyond 928
    n0 = 3*nhalf;
    n1 = 1:nmax;
    time_waveform(n1) = exp(-((n1-n0)/nhalf).^2);
    figure(1); plot(n1,time_waveform);
else       % sqaure wave
    imax = 200;
    nmax = 50;
end

Ez = zeros(imax,1); % Initializing the array
Hy = zeros(imax-1,1); % Initializing the array
%% TFSF parameters

TFSF = 100;
Hy_inc = zeros(imax-1,1);
Ez_inc = zeros(imax,1);


%%
for N = 1:nmax
%    

    for i = 1:imax-1  % regular update Hy_inc and Hy field
        Hy_inc(i) = Da*Hy_inc(i) + Db*(Ez_inc(i+1)-Ez_inc(i));
        Hy(i) = Da*Hy(i) + Db*(Ez(i+1)-Ez(i));
        
    end
    
     % correctons for Hy field for TFSF interfaces
     Hy(TFSF-1) = Da*Hy(TFSF-1) - Db*(Ez_inc(TFSF));
     Hy(imax-TFSF+1) = Da*Hy(imax-TFSF+1) + Db*Ez_inc(imax-TFSF+1);
        
     Ez_inc(1) = time_waveform(N);
     
    
    for i = 2:imax-1  % update Ez_inc and Ez field
        Ez_inc(i) = Ca*Ez_inc(i) + Cb*(Hy_inc(i)-Hy_inc(i-1));
        Ez(i) = Ca*Ez(i) + Cb*(Hy(i)-Hy(i-1));
    end
    
    %corrections for Ez field for TFSF interfaces
    Ez(TFSF) = Ca*Ez(TFSF) - Cb*Hy_inc(TFSF-1);
    Ez(imax-TFSF+1) = Ca*Ez(imax-TFSF+1) + Cb*Hy_inc(imax-TFSF+1);
    
    Ez(imax/2) = 0; % PEC where the crane recieving antenna is
    
    
    % plotting Ez
    figure(2);
    plot(Ez_inc,'LineWidth',2);
    set(gca, 'FontSize',14);
    axis([0 imax -1 1]);
    pause(0.001);
    title('gaussian pulse, S = 0.99, PEC center', 'FontSize', 14);
    ylabel('Ez(i)', 'FontSize', 14);
    xlabel('Grid i coordinate','FontSize', 14);
        
            
        
 end
    
    


