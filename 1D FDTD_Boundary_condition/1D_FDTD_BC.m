% this code is to test 1D FDTD simulation with a square pulse and a
% guassian source pulse with boundary conditions defined 
% starting from left |PEC -> absorbing material -> air  
% the interface between absorbing material and air is the PML layer
% To avoid reflections from the inciddent wave on the Simulated E field, we
% need a material with has zero reflection coefficient at the start( left
% side) and which can absorb. 
clc;
clear;
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
% common variables used in FDTD code 
eps = 8.85418782*1e-12; % m-3 kg-1 s4 A2 permittivity of free space
mu = 1.25663706*1e-6; % m kg s-2 A-2 permeability of free space
delta = 4.3*1e-3; % dx( for the avalanche case with ISM band 820Mhz - 980Mhz, with 10 grid/wavelength, wavelegnth = c/f_highest)
c = 2.99792458*1e+08;% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
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

% NEW VARIABLES FOR BOUNDARY CONDITION
% we introduce a polynomial grading for sigma(conductivity) in the
% absorbtion layer to avoid numerical reflections due to discretization
% errors
PML = 10; % thickness of the absrobing layer ( grid/cells)
m = 3; % order of the polynomial grading of sigma
eta = sqrt(mu/eps); % charecterisitc impedance of air
mu_r = 1; % air
eps_r = 1; % air
sigma_max = (0.8*(m+1))/(eta*delta*sqrt(mu_r*eps_r));

% Defining arrays for boundary conditions

Ez = zeros(imax,1); % Initializing the array
Hy = zeros(imax,1); % Initializing the array

sigma = zeros(imax,1);  % initializing array of sigma same size as Ez
sigma_star = zeros(imax,1); % initializing array of sigma same size as Ez

Ca = zeros(imax,1);
Cb = zeros(imax,1);
Da = zeros(imax,1);
Db = zeros(imax,1);


for i = 2:PML+1
    sigma(i) = (((PML+1.5-i)/(PML+0.5))^m)*sigma_max;
end

for i = 1:PML
    sigma_star(i) = (((PML+1-i)/(PML+0.5))^m)*sigma_max*(mu/eps);
end

% coefficients

for i = 1:imax
    Ca(i) = (1-((sigma(i)*dt)/(2*eps)))/(1+((sigma(i)*dt)/(2*eps)));
    Cb(i) = (dt/(eps*delta))/(1+((sigma(i)*dt)/(2*eps)));
    Da(i) = (1-((sigma_star(i)*dt)/(2*mu)))/(1+((sigma_star(i)*dt)/(2*mu)));
    Db(i) = (dt/(mu*delta))/(1+((sigma_star(i)*dt)/(2*mu)));   
end


% Da = 1; % constant from faraday's law
% Ca = 1; % constant from ampere's law
% Db = dt/(mu*delta); % constant from faraday's law
% Cb = dt/(eps*delta); % constant from ampere's law




%%
for N = 1:nmax
    if S1 == 0
        if N <= 20
            Ez(imax/4) = 1;
        else
            Ez(imax/4) = 0;
        end
    else
        
        Ez(imax/4) = sin(2*pi*f0*(N*dt-t0))*exp(-((N*dt-t0)/(thalf))^2);
    end
    
    
    for i = 1:imax-1  % update H field
        Hy(i) = Da(i)*Hy(i) + Db(i)*(Ez(i+1)-Ez(i));
        
    end
   % Ez(imax) = 0;
    for i = 1:imax-1  % update E field
        
        Ez(i+1) = Ca(i+1)*Ez(i+1) + Cb(i+1)*(Hy(i+1)-Hy(i));  % Ca and Cb are valid from i = 2 as sigma is not defined for i = 1
       
    end
    
    Ez_obs1(N) = Ez(imax/4 - 10);
    Ez_obs2(N) = Ez(imax/4 + 10); % reference observation point
    
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

