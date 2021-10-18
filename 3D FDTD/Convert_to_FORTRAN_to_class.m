%************************************************************************
% 1-D FDTD Code of a Gaussian Pulse Propagating in Free Space
%************************************************************************

% This is a finite-difference time-domain (FDTD) code that solves 
% Maxwell's equations in 1-D.  A plane wave propagates in free space
% along the x-direction.  The unknowns are are Ez and Hy.
% The source is a Gaussian waveform at the left-most position of the
% the grid (i = 1). 

clear all;

%*******************************
% DEFINE VARIABLES AND PARAMETERS
%*******************************

% Speed of light (m/s)
c = 2.99792456E8;

% permittivity of free space
eps = 8.85418782E-12;

% permeability of free space
mu = 4*pi*1E-7;

% space increment
dx = 1E-3;

% time step increment (S = 1)
S = 0.99;
dt = S*dx/c;

% imax is the number of grid cells in the x-direction
imax = 400;

% nmax is the total number of time steps to be run.
nmax = 410; 
                                     
% nhalf is the half-width of the Gaussian Pulse
nhalf = 20;
% no is how many time steps the center of the Gaussian is delayed
% in time from the start of the simulation
no = nhalf*3;

% arrays to hold the Ez and Hy fields
Ez = zeros(1,imax);
% start and end the grid on an Ez, so one less Hy value compared to Ez
Hy = zeros(1,imax-1); 

% Pre-calculate the updating coefficients 
Caz = 1.0;
Cbz = dt/(eps*dx);
Day = 1.0;
Dby = dt/(mu*dx);

%**************************************************************************
% Start the time-stepping loop
%**************************************************************************

for n=1:nmax;
   
%********************************
% update the Hy fields
%********************************
    for i = 1:imax-1
        Hy(i) = Day*Hy(i) + Dby*(Ez(i+1) - Ez(i));
    end
    
    
%********************************    
% Update the Ez fields 
%********************************

    % implement the source on Ez at i = 1
    Ez(1) = exp(-((n-no)/nhalf)^2);
  
    % Update the Ez fields 
    for i = 2:imax-1
        Ez(i)= Caz*Ez(i) + Cbz*(Hy(i) - Hy(i-1));
    end
 
%********************************    
% Output
%********************************

    % Plot the Ez values at each time step and label the axes.
    figure(1); plot(Ez,'LineStyle','-','LineWidth',2);
    set(gca,'FontSize',14);
    axis([0 400 -1.5 1.1]);
    title('Gaussian Pulse, Free Space', 'FontSize',16);
    xlabel('Grid i coordinate', 'FontSize',14);
    ylabel('Ez(i)', 'FontSize',14);
    pause(0.01)

end
%**************************************************************************
% End the time-stepping loop
%**************************************************************************

