% this code is to test 1D FDTD simulation with a square pulse and a
% guassian source pulse
clc;
clear;
close all;
%% Input S (courant stability factor)
n = input('Enter the parameter S(courant stability factor): ');
% ensures that numerical omega is real ( refer numerical dispersion
% relationship)

switch n
    case 1
        S = 0.5;
    case 2
        S = 0.99*(1/sqrt(2)); %2d courant stability 
    case 3
        S = 0.99;
        
end

%% Initializing variables
eps = 8.85418782*1e-12; % m-3 kg-1 s4 A2 permittivity of free space
mu = 1.25663706*1e-6; % m kg s-2 A-2 permeability of free space
c = 2.99792458*1e+08; % m/s speed of light
Da = 1; % constant from faraday's law
Ca = 1; % constant from ampere's law
imax = 100;
kmax = 100;
nmax = 100;
delta = 15;% in m => dx(gps abckup system antenna of frequency 1Mhz, with 20 grid/wavelength, wavelegnth = c/f_highest)% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
n1 = 1:nmax;
freq_antenna = 1*1e6;
Time_waveform(n1) = sin(((2*pi*(freq_antenna))*(n1)*dt));
figure(1); plot(n1,Time_waveform);
Db = dt/(mu*delta); % constant from faraday's law
Cb = dt/(eps*delta); % constant from ampere's law

Ez = zeros(imax,kmax); % Initializing the array
Hy = zeros(imax-1,kmax-1); % Initializing the array
Ex = zeros(imax,kmax); % Intializing the array
Hy_vert = zeros(nmax,1);

%%
for N = 1:nmax
   
    for k = 1:kmax-1
        for i = 1:imax-1  % update Hy field
            Hy(i,k) = Da*Hy(i,k) + Db*(Ez(i+1,k) - Ez(i,k) + Ex(i,k) - Ex(i,k+1));
            
        end
    end
    
    % adding source
    Hy(imax/2,kmax/2) = sin(((2*pi*(freq_antenna))*(N)*dt));
    
    for k = 1:kmax-1
        for i = 2:imax-1  % update Ez field
            Ez(i,k) = Ca*Ez(i,k) + Cb*(Hy(i,k) - Hy(i-1,k));
            
        end
    end
    
    for k = 2:kmax-1
        for i = 1:imax-1  % update Ex field
            Ex(i,k) = Ca*Ex(i,k) + Cb*(Hy(i,k-1)- Hy(i,k));
        end
    end
   
    % plotting Hy
   
    figure(2);
    pcolor(Hy');
    shading flat;
    axis('equal','tight');
    colorbar;
    set(gca, 'FontSize',14);
    pause(0.001);
    title('Free space grid', 'FontSize', 14);
    ylabel('Grid(i) coordinate', 'FontSize', 14);
    xlabel('Grid(k) coordinate','FontSize', 14);
   
    Hy_vert(N) = Hy(imax/2,kmax/2);
  
    
end


% figure(3);
% plot(Hy_vert);

