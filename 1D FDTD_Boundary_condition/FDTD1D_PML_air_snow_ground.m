% this code is to test 1D FDTD simulation with a square pulse and a
% guassian source pulse with boundary conditions defined
% starting from left |PEC -> absorbing material -> air -> snow -> human
% body -> snow -> ground
% the interface between absorbing material and air is the PML layer
% To avoid reflections from the inciddent wave on the Simulated E field, we
% need a material with has zero reflection coefficient at the start( left
% side) and which can absorb.
clc;
clear;
close all;
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
eps0 = 8.85418782*1e-12; % m-3 kg-1 s4 A2 permittivity of free space
mu0 = 1.25663706*1e-6; % m kg s-2 A-2 permeability of free space
delta = 4.3*1e-3; % dx( for the avalanche case with ISM band 820Mhz - 980Mhz, with 10 grid/wavelength, wavelegnth = c/f_highest)
c = 2.99792458*1e+08;% m/s speed of light
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
S1 = 1; % flag to choose square wave vs gaussian wave
if S1 ==1  % gaussian wave
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
eta = sqrt(mu0/eps0); % charecterisitc impedance of air
mu_r_air = 1; % air
eps_r_air = 1; % air
eps_r_snow = 1.6;
eps_r_muscle = 55.03;
eps_r_ground = 7;
sigma_r_ground = 0.01;
sigma_r_muscle = 0.99429;
sigma_r_snow = 0.00001;
sigma_max = (0.8*(m+1))/(eta*delta*sqrt(mu_r_air*eps_r_air));


% NEW VARIABLES to account for snow, person and ground interface
i_source = 30; % i = 30
i_surf = round(30 + 30/delta); % i = 7007
i_ground = round(i_surf + 20/delta);
imax = i_ground + 1300; % 8 skin depth deep ( half one way 4 skindepth = 4*1.4 = 5.6 => 5.6/delta = 1300)
nmax = 19400;

% Defining arrays for boundary conditions
Ez = zeros(imax,1); % Initializing the array
Hy = zeros(imax,1); % Initializing the array

sigma = zeros(imax,1);  % initializing array of sigma same size as Ez
sigma_star = zeros(imax,1); % initializing array of sigma same size as Ez

eps = zeros(imax,1);


Ca = zeros(imax,1);
Cb = zeros(imax,1);
Da = zeros(imax,1);
Db = zeros(imax,1);

% changing relative permeability for different interfaces ( sigma varies
% with distance in PML)

for i = 2:PML+1
    sigma(i) = (((PML+1.5-i)/(PML+0.5))^m)*sigma_max;
end

for i = 1:PML
    sigma_star(i) = (((PML+1-i)/(PML+0.5))^m)*sigma_max*(mu0/eps0);
end

body_flag = 1;
if body_flag == 1
    eps(1:i_surf) = eps_r_air*eps0;
    eps(i_surf+1:round(i_surf+(10/delta))) = eps_r_snow*eps0;
    eps(round(i_surf+(10/delta)+1):round(i_surf+(10.45)/(delta))) = eps_r_muscle*eps0;
    eps(round(i_surf+(10.45)/(delta)+1):i_ground) = eps_r_snow*eps0;
    eps(i_ground+1:imax) = eps_r_ground*eps0;
    
    sigma(i_surf+1:round(i_surf+(10/delta))) = sigma_r_snow;
    sigma(round(i_surf+(10/delta)+1):round(i_surf+(10.45)/(delta))) = sigma_r_muscle;
    sigma(round(i_surf+(10.45)/(delta)+1):i_ground) = sigma_r_snow;
    sigma(i_ground+1:imax) = sigma_r_ground;
    
else
    eps(1:i_surf) = eps_r_air*eps0;
    eps(i_surf+1:i_ground) = eps_r_snow*eps0;
    eps(i_ground+1:imax) = eps_r_ground*eps0;
    sigma(i_surf+1:i_ground) = sigma_r_snow;
    sigma(i_ground+1:imax) = sigma_r_ground;
    
end

figure(1);
plot((1:imax),eps,'LineWidth', 2);
title('Relative Permitivity vs grid points');
xlabel('grid(i)','LineWidth', 14);
ylabel('relative permitivitty', 'LineWidth', 14);

figure(2);
plot((1:imax),sigma,'LineWidth', 2);
title('Relative Permitivity vs grid points');
xlabel('grid(i)','LineWidth', 14)
ylabel('conductivity sigma (S/m)','LineWidth', 14);


% coefficients
%(Ca and Cb defined for point i is not relevant and is ignored in Ez calculation)
for i = 1:imax
    
    Ca(i) = (1-((sigma(i)*dt)/(2*eps(i))))/(1+((sigma(i)*dt)/(2*eps(i))));
    Cb(i) = (dt/(eps(i)*delta))/(1+((sigma(i)*dt)/(2*eps(i))));
    Da(i) = (1-((sigma_star(i)*dt)/(2*mu0)))/(1+((sigma_star(i)*dt)/(2*mu0)));
    Db(i) = (dt/(mu0*delta))/(1+((sigma_star(i)*dt)/(2*mu0)));
end



Ez_obs1 = zeros(nmax,1); % for finding relative error ( obs point 10 grid cells before source)
Ez_obs2 = zeros(nmax,1); % for finding relative error ( obs point 10 grid cells after source)
n1=1:nmax;
time_waveform = sin(2*pi*f0*(n1.*dt-t0)).*exp(-((n1.*dt-t0)/(thalf)).^2);



%%
for N = 1:nmax
    
    for i = 1:imax-1  % update H field
        Hy(i) = Da(i)*Hy(i) + Db(i)*(Ez(i+1)-Ez(i));
        
    end
    
    for i = 1:imax-1  % update E field
        
        Ez(i+1) = Ca(i+1)*Ez(i+1) + Cb(i+1)*(Hy(i+1)-Hy(i));  % Ca and Cb are valid from i = 2 as sigma is not defined for i = 1
        
    end
    
    if S1 == 0  %% adding source
        if N <= 20   % square wave
            Ez(imax/4) = 1;
        else
            Ez(imax/4) = 0;
        end
    else
        Ez(i_source) = Ez(i_source) - (dt/eps0).*time_waveform(N); % soft source
        
    end
    
    Ez_obs1(N) = Ez(i_source - 10);
    Ez_obs2(N) = Ez(i_source + 10); % reference observation point
    
    
    
    % plotting Ez
    if S1 == 0
        figure(3);
        plot(Ez,'LineWidth',2);
        set(gca, 'FontSize',14);
        title('Free space grid', 'FontSize', 14);
        ylabel('Ez(i)', 'FontSize', 14);
        xlabel('Grid i coordinate','FontSize', 14);
    else
        if mod(N,200) == 0
            figure(3);
            plot(Ez,'LineWidth',2);
            set(gca, 'FontSize',14);
            axis([0 imax -1 1]);
            %pause(0.001);
            title('Free space grid', 'FontSize', 14);
            ylabel('Ez(i)', 'FontSize', 14);
            xlabel('Grid i coordinate','FontSize', 14);
            
            
        end
    end
    
    
end
% define ylimits so we know how long to draw the lines
ylimits = ylim;
hold on;
% add line for snow surface
plot([i_surf,i_surf], [ylimits(1,1),ylimits(1,2)], '--k');
% add line for ground surface
plot([i_surf+i_ground,i_surf+i_ground], ...
    [ylimits(1,1),ylimits(1,2)], '--g');
% add line for body surface
plot([i_surf+round(10/delta),i_surf+round(10/delta)], ...
    [ylimits(1,1),ylimits(1,2)], '--r');
% add line for lower body surface
plot([i_surf+round((10+0.45)/delta),...
    i_surf+round((10+0.45)/delta)], ...
    [ylimits(1,1),ylimits(1,2)], '--r');
% add label for air region
text(3000,ylimits(1,2)*7/8, {'Air'}, 'FontSize',14);
% add label for snow region
text(7100,ylimits(1,2)*7/8, {'Snow'}, 'FontSize',14);
% add label for ground region
text(11550,ylimits(1,2)*7/8, {'Ground'},'Color','green','FontSize',14);
% add label for body region
text(8250,ylimits(1,2)*6/8, {'Body'}, 'Color','red', 'FontSize',14);
hold off;


figure(4);
plot(n1*dt*1e6,Ez_obs2,'LineWidth', 2);
title('Ez at observation point', 'LineWidth', 14);
xlabel('Time (microseconds)', 'LineWidth', 14);
ylabel('Ez', 'LineWidth', 14);

if body_flag == 1
    save Ez_obs_body.dat Ez_obs2 '-ascii';
    figure(5);
    plot(i_surf+round(10/delta):i_surf+round((10+0.45)/delta),Ez(i_surf+round(10/delta):i_surf+round((10+0.45)/delta)),'LineWidth', 2);
    title('Ez on body', 'LineWidth', 14);
    xlabel('grid(i)', 'LineWidth', 14);
    ylabel('Ez(i)', 'LineWidth', 14);
else
    save Ez_obs_nobody.dat Ez_obs2 '-ascii';
end
relative_error = abs(Ez_obs1-Ez_obs2)./abs(max(Ez_obs2));


figure(6);
plot(Ez(i_ground:imax),'LineWidth', 2);
title('Ez at ground showing decay', 'LineWidth', 14);
xlabel('grid(i)', 'LineWidth', 14);
ylabel('Ez(i)', 'LineWidth', 14);





