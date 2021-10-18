% this code is to test 2D FDTD simulation with a sinewave pulse and a
% assuming reflections between ground and inosphere waveguide and PML. Also
% add SIBC on the ground surface

clc;
clear;
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
reflection_no = 1;
if reflection_no == 1
    imax = 750;
    %PML = 100; % thickness of the absrobing layer ( grid/cells)
    
else
    imax = 160;
    %PML = 10;
end

% imax = 160;
nmax = 2000;
delta = 650;% in m =>ength, the tallest antenna possible
kmax = round(78000/delta); % ( altitude of reflection height from ionosphere and to find grid cells divide by delta)
dt = (S*delta)/c; % dt defined from lecture, as dt -> 0, the error goes down
n1 = 1:nmax;
freq_antenna = 1.0*1e4; % for wave to reflect the antenna freq should be in VLF range (3-30Khz)
Time_waveform(n1) = sin(((2*pi*(freq_antenna))*(n1)*dt));
figure(1); plot(n1,Time_waveform);


%% Ionosphere profile
sigma = zeros(kmax-1,1);
electron_charge = -1.6E-19;  % electron charge
electron_mass = 9.1E-31;     % electron mass
sharpness = 0.5;             % sharpness of the ionosphere profile
reflect_altitude = 75.0;     % reflection height of the ionosphere day time (75km and night = 85km)
Ca = zeros(kmax-1,1);
Cb = zeros(kmax-1,1);
Ca_1 = zeros(kmax-1,1);
Cb_1 = zeros(kmax-1,1);

for k = 1:kmax-1
    collision_freq = 1.816e11 * exp(-0.15 * k*delta/1000);
    electron_density = 1.43E13 * exp(-0.15 * reflect_altitude) * ...
        exp((sharpness-0.15)*(k*delta/1000 - reflect_altitude));
    sigma(k) = (electron_density * electron_charge^2) / ...
        (electron_mass * collision_freq);
    Ca(k) = (1-((sigma(k)*dt)/(2*eps0)))/(1+((sigma(k)*dt)/(2*eps0)));
    Cb(k) = (dt/(eps0*delta))/(1+((sigma(k)*dt)/(2*eps0)));
    %     Ca_1(k) = (1-((dt*sigma(k))/(2*eps0)))/(1+((dt*sigma(k))/(2*eps0))); % constant from ampere's law
    %     Cb_1(k) = ((4*dt)/(eps0*delta))/(1+((dt*sigma(k))/(2*eps0)));
end


%% Adding the PML layer to the right along X direction

% NEW VARIABLES FOR BOUNDARY CONDITION
% we introduce a polynomial grading for sigma(conductivity) in the
% absorbtion layer to avoid numerical reflections due to discretization
% errors
PML = 100;
m = 3; % order of the polynomial grading of sigma
eta = sqrt(mu0/eps0); % charecterisitc impedance of air
mu_r = 1; % air
eps_r = 1; % air
sigma_max = (0.8*(m+1))/(eta*delta*sqrt(mu_r*eps_r));

sigma_E_PML = zeros(PML,1);
b_E = zeros(PML,1);
c_E = zeros(PML,1);
Ca_PML = zeros(PML,1);
Cb_PML = zeros(PML,1);
psi_Ezx = zeros(PML+1,kmax-1);

sigma_H_PML = zeros(PML-1,1);
b_H = zeros(PML-1,1);
c_H = zeros(PML-1,1);
psi_Hyx = zeros(PML,kmax-1);


for i = 2:PML+1
    sigma_E_PML(i) = (((PML+1.5-i)/(PML+0.5))^m)*sigma_max;
    % sigma_E_PML(i) = sigma_max*((PML-(i-1.5))/(PML-1.0))^m;
    b_E(i) = exp(-((sigma_E_PML(i)))*dt/eps0);
    c_E(i) = exp(-((sigma_E_PML(i)))*dt/eps0) - 1;
end

for i = 1:PML
    sigma_H_PML(i) = (((PML+1-i)/(PML+0.5))^m)*sigma_max;
    %sigma_H_PML(i) = (((PML+1-i)/(PML-1))^m)*sigma_max;
    b_H(i) = exp(-((sigma_H_PML(i)))*dt/eps0);
    c_H(i) = exp(-((sigma_H_PML(i)))*dt/eps0) - 1;
end




%% SIBC layer
sigma_ground = 0.01;
sigma_ocean = 0.33;
sigma_air = 0;
eps_ground = eps0*1.5;

Rs_f = sqrt((mu0*2*pi*freq_antenna)/(2*sigma_ground)); % using the surface impedance bounary condition ( dissipation with Rs)
Ls_f = sqrt((mu0)/(2*2*pi*freq_antenna*sigma_ground));  % using the surface impedance bounary condition ( inductance with Ls)
DA_SIBC = ((mu0*delta) - (0.5*dt*Rs_f) + Ls_f)/((mu0*delta) + (0.5*dt*Rs_f) + Ls_f); % constant for updating the Hy field for ground propagation
DB_SIBC = (dt)/((mu0*delta) + (0.5*dt*Rs_f) + Ls_f);

%%

Da = 1; % constant from faraday's law
Db = dt/(mu0*delta); % constant from faraday's law

Ez = zeros(imax,kmax-1); % Initializing the array
Hy = zeros(imax-1,kmax-1); % Initializing the array
Ex = zeros(imax-1,kmax); % Intializing the array
Ezmax = zeros(imax,kmax-1);
period = round((1/freq_antenna)*(1/dt));
Hy_vert = zeros(nmax,1);
Ez_obs = zeros(nmax,1);
%%
for N = 1:nmax
    
    for k = 1:kmax-1
        if k == 1
            for i = 1:imax-1
                Hy(i,k) = DA_SIBC*Hy(i,k) + DB_SIBC*(Ez(i+1,k) - Ez(i,k) - Ex(i,k+1)); % adding the effect of ground propagation
            end
        else
            for i = 1:imax-1  % regular update Hy field
                Hy(i,k) = Da*Hy(i,k) + Db*(Ez(i+1,k) - Ez(i,k) + Ex(i,k) - Ex(i,k+1));
                
            end
        end
        
        
    end
    
    % right side PML layer psi_Hyx update
    
    for k = 1:kmax-1
        for i = imax-PML:imax-1
            psi_Hyx(imax-i,k) = b_H(imax-i)*psi_Hyx(imax-i,k) + c_H(imax-i)*(Ez(i+1,k) - Ez(i,k))*(1/delta);
            Hy(i,k) = Da*Hy(i,k) + dt/mu0*psi_Hyx(imax-i,k);
        end
        
    end
    
    for k = 2:kmax-1
        for i = 1:imax-1  % update Ex field
            Ex(i,k) = Ca(k)*Ex(i,k) + Cb(k)*(Hy(i,k-1)- Hy(i,k));
        end
    end
    
    
    for k = 1:kmax-1
        for i = 2:imax-1  % update Ez field
            Ez(i,k) = Ca(k)*Ez(i,k) + Cb(k)*(Hy(i,k) - Hy(i-1,k));
            
        end
    end
    
    % right side PML layer for Ez propagation psi_Ezx
    for k = 1:kmax-1
        for i = imax-PML:imax-1
            psi_Ezx(imax-i+1,k) = b_E(imax-i+1)*psi_Ezx(imax-i+1,k) + c_E(imax-i+1)*(Hy(i,k) - Hy(i-1,k))/delta;
            %Ez(i,k) = Ca(k)*Ez(i,k) + Cb(k)*psi_Ezx(imax-i+1,k);
            Ez(i,k) = Ez(i,k) + dt/eps0*psi_Ezx(imax-i+1,k);
        end
    end
    
    
    
    for k = 2:kmax-1     % update Ez for source symmetric along left side of the grid using ampere' law in integral form
        Ca_1(k) = (1-((dt*sigma(k))/(2*eps0)))/(1+((dt*sigma(k))/(2*eps0)));
        Cb_1(k) = ((4*dt)/(eps0*delta))/(1+((dt*sigma(k))/(2*eps0)));
        Ez(1,k) = Ca_1(k)*Ez(1,k) + Cb_1(k)*Hy(1,k);
    end
    
    
    
    % adding source
    Ez(1,1) = sin(((2*pi*(freq_antenna))*(N)*dt));
    
    % loop for adding Ezmax
    for k = 1:kmax-1
        for i = 1:imax-1
            if (Ezmax(i,k) < abs(Ez(i,k)) && N > (nmax -2*period))
                Ezmax(i,k) = abs(Ez(i,k));
            end
        end
        
    end
    
    
    Ez_obs(N) = Ez(140,1);
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
        
        
        
        
    end
end

% plotting Ezamx on logscale
figure(3);
subplot(2,1,2);
pcolor(log10(abs(Ezmax')));
caxis([-3.2 -2.0])
shading flat;
axis('equal','tight');
colorbar;
set(gca, 'FontSize',14);
title('Free space grid with inosphere and ground along vertical axis', 'FontSize', 14);
ylabel('Grid(i) coordinate', 'FontSize', 14);
xlabel('Grid(k) coordinate','FontSize', 14);
hold on;

% looking closely at Ezmax along surface of earth
figure(4); plot(log10(Ezmax(1:imax-PML,1)), 'LineWidth',2);
xlabel('Grid i coordinate');
ylabel('Log scale plot of Max Ez')
title('maximum Ez at steady state w.r.t earth surface')
hold on;


%% relative error

if reflection_no == 1
    save Ez_obs_large.dat Ez_obs '-ascii';
    figure(4);
    plot(n1,Ez_obs,'LineWidth', 2);
    title('Ez at observation point', 'LineWidth', 14);
    xlabel('Time (microseconds)', 'LineWidth', 14);
    ylabel('Ez', 'LineWidth', 14);
else
    save Ez_obs_small.dat Ez_obs '-ascii';
    figure(4);
    plot(n1,Ez_obs,'LineWidth', 2);
    title('Ez at observation point', 'LineWidth', 14);
    xlabel('Time (microseconds)', 'LineWidth', 14);
    ylabel('Ez', 'LineWidth', 14);
end




load Ez_obs_large.dat
load Ez_obs_small.dat

relative_error = abs((Ez_obs_large)-(Ez_obs_small))./abs(max(Ez_obs_small));
figure(5);
semilogy(relative_error);

% figure(3);
% plot(Hy_vert);

