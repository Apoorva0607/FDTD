%*******************************
% Homework 3 Example Solution
%*******************************

% This is a finite-difference time-domain (FDTD) code that solves 
% Maxwell's equations in 1-D.  A plane wave propagates 
% in the x-direction.  The unknowns are are Ez and Hy.
% The source is a Gaussian waveform at the left-most position of the
% the grid (i = 1). A perfect absorbing boundary condition (ABC) is 
% implemented on Ez at i = imax (assuming S = 1).

clear all;
close all;
format long;

%*******************************
% DEFINE VARIABLES AND PARAMETERS
%*******************************

% Speed of light (m/s)
c = 2.99792456E8;

% permittivity of free space
epso = 8.85418782E-12;

% permeability of free space
muo = 4*pi*1E-7;

% imax is the number of grid cells in the x-direction
imax = 750;  % long enough to see big null
kmax = 120; %model only up to 76800 km or else wavelength too small and PML goes unstable
% wavelength is 755 m, skin depth is 1.55 km at kmax = 128 (3 cells would
% be better, but getting at least two.


% Nmax is the total number of time steps to be run.
nmax = 2000;  % we want to run it long enough to get steady state over all grid cells of interest 
%if we assume things propagate at the speed of light straight across the
%grid, it will take imax*delta / 3E8 / dt time steps to get to the end
% but propagating at an angle, so let's assume 
                         
freq = 1E4;
% nhalf is the half-width of the Gaussian Pulse
%nhalf = 100;
% no is how many time steps the center of the Gaussian is delayed
% in time from the start of the simulation
%no = nhalf*3;

% space increment
delta = 650;


% time step increment (S = 1)
S = 0.9;  % Had to lower to make cylindrical axis work.
dt = S*delta/c/sqrt(2);

 sigma = zeros(kmax-1,1);
 electron_charge =-1.0*1.6E-19;
electron_mass=9.1E-31;
      beta=0.5;
      h_init=75.0;


 for k = 1:kmax-1
     nu = 1.816e11 * exp(-0.15*k*delta/1000);  
     Ne = 1.43E13*exp(-0.15*h_init)*exp((beta-0.15)*(k*delta/1000-h_init));
     sigma(k)=(Ne * electron_charge^2) / (electron_mass * nu);
 end
 
% sigma_star = zeros(imax-1,1);

%**********************************************
% PML 
%**********************************************

%thickness of the PML in grid cells (Ez fields), including the last PEC Ez
nPML = 100;
% the sigma losses use a polynomial grading with m = 3 (see Equation 7.60)
m = 3;
ma = 1;
% sigma max is equal to sigma opt (Equation 7.67)
sigma_max = 0.8*(m+1)/(sqrt(muo/epso)*delta*sqrt(1.0));
alpha_max = 0.0;%0.24;
kappa_max = 1.0;

sigma_PML = zeros(nPML,1);
sigma_star_PML = zeros(nPML-1,1);
alphae = zeros(nPML,1);
alphah = zeros(nPML-1,1);
kappae = zeros(nPML,1);
kappah = zeros(nPML-1,1);
be = zeros(nPML,1);
bh = zeros(nPML-1,1);
ce = zeros(nPML,1);
ch = zeros(nPML-1,1);

psi_Ezx_2 = zeros(nPML+1,kmax-1);
psi_Hyx_2 = zeros(nPML,kmax-1);

%eta = sqrt(muo/epso/7)*(1-j*0.01/(2*pi*900E6*epso*7))^(-1/2);

%for i = imax-nPML:imax-1
    %sigma(i) = sigma_max*((i-(imax-nPML)+0.5)/(nPML-1.0))^m;  %%% add on ground sigma?
for i = 2:nPML+1  % can reuse all of these on all sides
    sigma_PML(i) = sigma_max*((nPML-(i-1.5))/(nPML-1.0))^m;  % same as for 1-D code on left side;  will reuse in reverse order on right side; 
    % call sigma_PML to distinguish from PML in main grid
    alphae(i) = alpha_max*((i-1.5)/(nPML-1.0))^ma;
    kappae(i) = 1.0+(kappa_max-1.0)* ((nPML-(i-1.5)) / (nPML - 1.0))^m;
    be(i) = exp(-(sigma_PML(i) / kappae(i) + alphae(i))*dt/epso);
    ce(i) = sigma_PML(i)*(be(i)-1.0)/(sigma_PML(i)+kappae(i)*alphae(i))/ kappae(i);
end
% Assign the sigma_star at the Hy and Hx locations (see Equation 7.60) 
% and Equation 7.8
% PML on Hy's goes from i = imax-nPML to i = imax-1
%for i = imax-nPML:imax-1
    %sigma_star(i) = sigma_max*((i-(imax-nPML)+1)/(nPML-1.0))^m *muo/epso;
%end
for i = 1:nPML % can reuse on all sides
    sigma_star_PML(i) = (sigma_max*((nPML-(i-1.0))/(nPML-1.0))^m);  % no multiply by mu/eps!!!!!!
    alphah(i) = alpha_max*((i-1.0)/(nPML-1.0))^ma; % max value at interface and decays to zero on edges p. 297
    kappah(i) = 1.0+(kappa_max-1.0)* ((nPML-(i-1.0)) / (nPML - 1.0))^m;
    bh(i) = exp(-(sigma_star_PML(i) / kappah(i) + alphah(i))*dt/epso);
    ch(i) = sigma_star_PML(i)*(bh(i)-1.0)/ (sigma_star_PML(i)+kappah(i)*alphah(i)) / kappah(i);
end

% arrays to hold the Ez and Hy fields
Ez = zeros(imax,kmax-1); %still have imax in x-direction as for 1-D code
% start and end the grid on an Ez, so one less Hy value compared to Ez
Ez_obs = zeros(nmax,1);
Ex = zeros(imax-1,kmax);
Ez_max = zeros(imax,kmax-1);
Hy = zeros(imax-1,kmax-1);
% Pre-calculate the updating coefficients (so don't have to recalculate 
% every time step). 
% Ca = 1;
% Cb = dt/(epso*delta);  % must have delta here or goes unstalbe == precision issue?
% Da = 1;
% Db = dt/(muo*delta);
Ca=zeros(kmax-1,1);
Cb=zeros(kmax-1,1);

% Ca(1:kmax-1) = 1.0;
% Cb(1:kmax-1) = dt / (epso*delta);
for k = 1:kmax-1
Ca(k) = (1-((sigma(k)*dt)/(2*epso)))./(1+((sigma(k)*dt)/(2*epso)));
Cb(k) = (dt/(epso*delta))/(1+((sigma(k)*dt)/(2*epso)));
end
Da = 1.0;
Db = dt/(muo*delta);

%Only need 1-D
%den_Ez = zeros(imax,1);
%den_Hy = zeros(imax-1,1);

%den_Ez(:) = 1.0;%delta;
%den_Hy(:) = 1.0;%delta;

% for i = imax-nPML:imax-1
%     denominator_Ez(i) = 1.0/(delta*kappae(i-(imax-nPML)+1));
% end
% for i = imax-nPML:imax-1
%     denominator_Hy(i) = 1.0/(delta*kappah(i-(imax-nPML)+1));
% end

%*******************************
% Start time-stepping loop
%*******************************

for n=1:nmax
    
        % Update the Hy fields (see Lecture 4d)
    for k = 1:kmax-1
       for i = 1:imax-1
          Hy(i,k) = Da*Hy(i,k) + Db*(Ez(i+1,k) - Ez(i,k) + Ex(i,k) - Ex(i,k+1));
       end
    end
          %Hy(1:imax-1,1:jmax-1) = Da*Hy(1:imax-1,1:jmax-1) + Db*(Ez(1:imax-1,1:jmax-1) - Ez([1:imax-1]+1,1:jmax-1)+Ex([1:imax-1],[1:jmax-1]+1) - Ex(1:imax-1,1:jmax-1));
         %PML
    for k = 1:kmax-1
        for i = imax-nPML:imax-1
            psi_Hyx_2(imax-i,k) = bh(imax-i)*psi_Hyx_2(imax-i,k) + ch(imax-i) *(Ez(i+1,k) - Ez(i,k))/delta;
            Hy(i,k) = Hy(i,k) + dt/muo*psi_Hyx_2(imax-i,k);
        end
    end
    
    
        % Update the Ex fields (see Lecture 4d)
      for k = 2:kmax-1
        for i = 1:imax-1
            Ex(i,k)= Ca(k)*Ex(i,k) + Cb(k)*(Hy(i,k-1) - Hy(i,k));
        end
       end
    
    % Update the Ez fields (see Lecture 4d)
    for k = 1:kmax-1
        for i = 2:imax-1
            Ez(i,k)= Ca(k)*Ez(i,k) + Cb(k)*(Hy(i,k) - Hy(i-1,k));
            %Ez(2:imax-1,1:jmax-1)= Ca*Ez(2:imax-1,1:jmax-1) + Cb*(Hy([2:imax-1]-1,1:jmax-1) - Hy(2:imax-1,1:jmax-1));
            if (Ez_max(i,k) < abs(Ez(i,k)) && n > (nmax-100))   % one period is 31.4 time steps (1/frequency /dt) so this ensures sampling over the last few periods.
                Ez_max(i,k) = abs(Ez(i,k));
            end
        end
    end
    

    for k = 2:kmax-1
Ca_azim = (1-((sigma(k)*dt)/(2*epso)))./(1+((sigma(k)*dt)/(2*epso)));
Cb_azim = (dt/(epso))/(1+((sigma(k)*dt)/(2*epso)));
        Ez(1,k) = Ca_azim*Ez(1,k)+ Cb_azim*2.0*pi*delta/2*(Hy(1,k))/(pi*(delta/2)^2);
    end
    
   % Ca_azim = 1.0;
%Cb_azim = dt / epso / delta;
 %       Ez(1,k) = Ca_azim*Ez(1,k)+ Cb_azim*(Hy(1,k));   % assume it is a rectangular cell!!!!!!
 %   end
    
    
        Ez(1,1) = sin(2*pi*n*freq*dt);
%     
    %PML
    for k = 1:kmax-1
        for i = imax-nPML:imax-1
            psi_Ezx_2(imax-i+1,k) = be(imax-i+1)*psi_Ezx_2(imax-i+1,k) + ce(imax-i+1) *(Hy(i,k) - Hy(i-1,k))/delta;
            Ez(i,k) = Ez(i,k) + dt/epso*psi_Ezx_2(imax-i+1,k);
        end
    end
    
    % Set the source 
    %Ez(i) = 1.0;
    %Ez(imax/4) = exp(-((n-no)/nhalf)^2)*sin(2*pi*(n-no)*freq*dt);
    
    Ez_obs(n) = Ez(140,1);
    
    if (mod(n,100) == 0) 
    % Plot the Ez values at each time step and label the axes.
    figure(1);
    pcolor((log10(abs(Ez')))); %, 'color','g','LineStyle','-','LineWidth',2);
    caxis([-3.2 -2.0])
    shading flat;
    axis('equal','tight');
    colorbar;
    set(gca, 'FontSize',14);
    pause(0.001);
   % axis([0 imax -1 1]);
   % title('Gaussian Pulse: S=1.0', 'FontSize',16);
    %xlabel('Grid i coordinate', 'FontSize',14);
    %ylabel('Ez(i)', 'FontSize',14);
 
        
    end

end
% End time-stepping

figure(11);
save Ez_obs_large.dat Ez_obs '-ascii';
figure(4);
n1 = 1:nmax;
plot(Ez_obs,'LineWidth', 2);
title('Ez at observation point', 'LineWidth', 14);
xlabel('Time (microseconds)', 'LineWidth', 14);
ylabel('Ez', 'LineWidth', 14);

figure(3);
pcolor(log10(abs(Ez_max'))); 
caxis([-3.2 -2.0])
shading flat;
axis('equal','tight');
colorbar;
set(gca, 'FontSize',14);
title('Free space grid with inosphere and ground along vertical axis', 'FontSize', 14);
ylabel('Grid(i) coordinate', 'FontSize', 14);
xlabel('Grid(k) coordinate','FontSize', 14);

%      pcolor((20.*log10(abs(Ez_max'./(max(max(Ez_max))))))); %, 'color','g','LineStyle','-','LineWidth',2);
%     shading flat;
%     colorbar;
%     set(gca,'FontSize',14);
%     caxis([-60 -20]);
%    % axis([0 imax -1 1]);
%    % title('Gaussian Pulse: S=1.0', 'FontSize',16);
%     %xlabel('Grid i coordinate', 'FontSize',14);
%     %ylabel('Ez(i)', 'FontSize',14);
    
%      figure(2); 
%      plot([1:imax-nPML].*delta,20.*log10(Ez_max(1:imax-nPML,1)./max(max(Ez_max)))) %don't plot end because will not
%      %see much of what is going on
%      set(gca,'FontSize',14);
%      title('Gaussian Pulse: S=1.0', 'FontSize',16);
%      xlabel('Grid i coordinate', 'FontSize',14);
%      ylabel('Ez(i)', 'FontSize',14);

