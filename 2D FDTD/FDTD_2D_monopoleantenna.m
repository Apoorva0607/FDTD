% this code is to test 2D FDTD simulation with a sinewave pulse and a
% assuming reflections between ground and inosphere.

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
        S = 0.9*(1/sqrt(2)); %2d courant stability for monopole antenna
    case 3
        S = 0.99;
        
end



%% Initializing variables
eps0 = 8.85418782*1e-12; % m-3 kg-1 s4 A2 permittivity of free space
mu0 = 1.25663706*1e-6; % m kg s-2 A-2 permeability of free space
c = 2.99792458*1e+08; % m/s speed of light
imax = 500;
nmax = 500;
delta = 650;% in m =>ength, the tallest antenna possible
kmax = 78000/delta; % ( altitude of reflection height from ionosphere and to find grid cells divide by delta)
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
n1 = 1:nmax;
freq_antenna = 1*1e4; % for wave to reflect the antenna freq should be in VLF range (3-30Khz)
Time_waveform(n1) = sin(((2*pi*(freq_antenna))*(n1)*dt));
figure(1); plot(n1,Time_waveform);


%% Ionosphere profile
sigma = zeros(kmax-1,1);
electron_charge = -1.6E-19;  % electron charge
electron_mass = 9.1E-31;     % electron mass
sharpness = 0.5;             % sharpness of the ionosphere profile
reflect_altitude = 75.0;     % reflection height of the ionosphere

Ca = zeros(imax,kmax);
Cb = zeros(imax,kmax);
Ca_1 = zeros(imax,kmax);
Cb_1 = zeros(imax,kmax);

for k = 1:kmax-1
    collision_freq = 1.816e11 * exp(-0.15 * k*delta/1000);
    electron_density = 1.43E13 * exp(-0.15 * reflect_altitude) * ...
        exp((sharpness-0.15)*(k*delta/1000 - reflect_altitude));
    sigma(k) = (electron_density * electron_charge^2) / ...
        (electron_mass * collision_freq);
    Ca(1,k) = (1-((sigma(k)*dt)/(2*eps0)))/(1+((sigma(k)*dt)/(2*eps0)));
    Cb(1,k) = (dt/(eps0*delta))/(1+((sigma(k)*dt)/(2*eps0)));
    Ca_1(1,k) = (1-((dt*sigma(k))/(2*eps0)))/(1+((dt*sigma(k))/(2*eps0))); % constant from ampere's law
    Cb_1(1,k) = ((4*dt)/(eps0*delta))/(1+((dt*sigma(k))/(2*eps0)));
end

Ca = repmat(Ca(1,:),[imax 1]);
Cb = repmat(Cb(1,:),[imax 1]);
Ca_1 = repmat(Ca_1(1,:),[imax 1]);
Cb_1 = repmat(Cb_1(1,:),[imax 1]);

%%
sigma_ground = 0.01;
sigma_air = 0;
eps_ground = eps0*1.5;

Da = 1; % constant from faraday's law
Db = dt/(mu0*delta); % constant from faraday's law

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
    
    
    for k = 1:kmax-1
        for i = 2:imax-1  % update Ez field
            Ez(i,k) = Ca(i,k)*Ez(i,k) + Cb(i,k)*(Hy(i,k) - Hy(i-1,k));
            
        end
    end
    
    for k = 1:kmax-1     % update Ez for source at first i=1,k=1 point using ampere' law in integral form
        Ez(1,k) = Ca_1(i,k)*Ez(1,k) + Cb_1(i,k)*Hy(1,k);
    end
    
    
    for k = 2:kmax-1
        for i = 1:imax-1  % update Ex field
            Ex(i,k) = Ca(i,k)*Ex(i,k) + Cb(i,k)*(Hy(i,k-1)- Hy(i,k));
        end
    end
    
    % adding source
    
    Ez(1,1) = sin(((2*pi*(freq_antenna))*(N)*dt));
    
    
    % plotting Ez
    if mod(N,100) == 0
        figure(2);
        pcolor(log10(abs(Ez')));
        caxis([-3.2 -2.0])
        shading flat;
        axis('equal','tight');
        colorbar;
        set(gca, 'FontSize',14);
        pause(0.001);
        title('Free space grid with inosphere and ground along vertical axis', 'FontSize', 14);
        ylabel('Grid(i) coordinate', 'FontSize', 14);
        xlabel('Grid(k) coordinate','FontSize', 14);
        
        Hy_vert(N) = Hy(imax/2,kmax/2);
        
    end
end


% figure(3);
% plot(Hy_vert);

