% this code is to test 1D FDTD simulation with a square pulse and a
% guassian source pulse
clc;
clear;
clear all;
%% Input S (courant stability factor)
n = input('Enter the parameter S(courant stabilityb factor): ');
% ensures that numerical omega is real ( refer numerical dispersion
% relationship)

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
delta = 4.3*1e-3; % dx( for the avalanche case with ISM band 820Mhz - 980Mhz, with 10 grid/wavelength, wavelegnth = c/f_highest)
c = 2.99792458*1e+08;% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
Db = dt/(mu*delta); % constant from faraday's law
Cb = dt/(eps*delta); % constant from ampere's law
S1 = 1; % flag to choose square wave vs gaussian wave
if S1 ==1  % gaussian wave
    imax = 25000; % need more grid point to compute field for gaussian wave
    nmax = 20000; % found by plotting the gaussian function with time
    f0 = 915*1e6; % since the ISM band is from 820MHz to 980MHz
    thalf = 2/(pi*13*1e6);% ieally should be 2/(pi*26MHz) but this would make the spectrum go beyond 928
    t0 = 3*thalf;
else       % sqaure wave
    imax = 200;
    nmax = 50;
end

Ez = zeros(imax,1); % Initializing the array
Hy = zeros(imax,1); % Initializing the array



%%
for N = 1:nmax
%     if S1 == 0
%         if N <= 20
%             Ez(imax/4) = 1;
%         else
%             Ez(imax/4) = 0;
%         end
%     else
%         
%         Ez(imax/4) = sin(2*pi*f0*(N*dt-t0))*exp(-((N*dt-t0)/(thalf))^2);
%     end
    
    
    for i = 1:imax-1  % update H field
        Hy(i) = Da*Hy(i) + Db*(Ez(i+1)-Ez(i));
        
    end
    Ez(imax) = 0;
    for i = 1:imax-1  % update E field
        
        Ez(i+1) = Ca*Ez(i+1) + Cb*(Hy(i+1)-Hy(i));
    end
    
    if S1 == 0
        if N <= 20
            Ez(imax/4) = 1;
        else
            Ez(imax/4) = 0;
        end
    else
        
        Ez(imax/4) = sin(2*pi*f0*(N*dt-t0))*exp(-((N*dt-t0)/(thalf))^2);
    end

    
    % plotting Ez
    if S1 == 0
        figure(2);
        plot(Ez,'LineWidth',2);
        set(gca, 'FontSize',14);
        title('Free space grid', 'FontSize', 14);
        ylabel('Ez(i)', 'FontSize', 14);
        xlabel('Grid i coordinate','FontSize', 14);
    else
        if mod(N,200) == 0
            figure(2);
            plot(Ez,'LineWidth',2);
            set(gca, 'FontSize',14);
            axis([0 imax -1 1]);
            pause(0.001);
            title('Free space grid', 'FontSize', 14);
            ylabel('Ez(i)', 'FontSize', 14);
            xlabel('Grid i coordinate','FontSize', 14);
        end
    end
    
    
end

