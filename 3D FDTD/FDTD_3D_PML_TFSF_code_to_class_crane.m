%************************************************************************
% 3-D FDTD Model for the Construction Site Design Challenge
%************************************************************************

% This is a finite-difference time-domain (FDTD) code that solves 
% Maxwell's equations in 3-D.  A plane wave propagates across the grid
% in the x-direction using the total-field scattered field formulation 
% from i = TFSF to i = imax-TFSF, and j = TFSF to j = jmax-TFSF and 
% from k = 1 to k = kmax-TFSF.  
% The grid is bounded by a perfectly matched layer (PML) that is "PML" 
% grid cells thick on the left, right, front, back, and top sides.

clear all;
%format long;

%**************************************************************************
% DEFINE VARIABLES AND PARAMETERS
%**************************************************************************

% Speed of light (m/s)
c = 2.99792456E8;

% permittivity of free space (F/m)
epso = 8.85418782E-12;

% permeability of free space (H/m)
muo = 4*pi*1E-7;

% imax is the number of grid cells in the x-direction, jmax is the number
% of grid cells in the y-direction, and kmax is the number of grid cells in
% the z-direction.
imax = 60;
jmax = 60;
kmax = 90; % 60 initially

% nmax is the total number of time steps 
nmax = 3000;  
          
% delta is the space increment in all three Cartesian directions (cubic
% cells)
delta = 1;

% The source is a sinusoid set to a center frequency of freq
%freq = c/(28*delta);  % 20 cellsl per wavelength
freq = 1280*1E3;

% time step increment (S = 1)
S = 0.99/sqrt(3);
dt = S*delta/c;

%**************************************************************************
% PML (Convolution PML with alpha = 0 and kappa = 1)
%**************************************************************************

% thickness of the PML in grid cells
PML = 10;
% the sigma losses use a polynomial grading with m = 3 
m = 3;
% sigma max is equal to sigma opt 
sigma_max = 0.8*(m+1)/(sqrt(muo/epso)*delta*sqrt(1.0));

sigma_E_PML = zeros(PML+1,1);
sigma_H_PML = zeros(PML,1);

bE = zeros(PML+1,1);
bH = zeros(PML,1);
cE = zeros(PML+1,1);
cH = zeros(PML,1);

% assign the sigmas in the PML region
for i = 2:PML+1  % can reuse all of these on all sides
    sigma_E_PML(i) = ((PML+1.5-i)/(PML+0.5))^m * sigma_max;
    bE(i) = exp(-sigma_E_PML(i)*dt/epso);
    cE(i) = bE(i)-1.0;
end

for i = 1:PML 
    sigma_H_PML(i) = ((PML+1-i)/(PML+0.5))^m * sigma_max;
    bH(i) = exp(-sigma_H_PML(i)*dt/epso);
    cH(i) = bH(i)-1.0;
end

% x-direction PML
psi_Ezx_1 = zeros(PML+1,jmax,kmax-1);
psi_Hyx_1 = zeros(PML,jmax,kmax-1);

psi_Eyx_1 = zeros(PML+1,jmax-1,kmax);
psi_Hzx_1 = zeros(PML,jmax-1,kmax);

psi_Ezx_2 = zeros(PML+1,jmax,kmax-1);
psi_Hyx_2 = zeros(PML,jmax,kmax-1);

psi_Eyx_2 = zeros(PML+1,jmax-1,kmax);
psi_Hzx_2 = zeros(PML,jmax-1,kmax);

% y-direction PML
psi_Exy_1 = zeros(imax-1,PML+1,kmax);
psi_Hzy_1 = zeros(imax-1,PML,kmax);

psi_Ezy_1 = zeros(imax,PML+1,kmax-1);
psi_Hxy_1 = zeros(imax,PML,kmax-1);

psi_Exy_2 = zeros(imax-1,PML+1,kmax);
psi_Hzy_2 = zeros(imax-1,PML,kmax);

psi_Ezy_2 = zeros(imax,PML+1,kmax-1);
psi_Hxy_2 = zeros(imax,PML,kmax-1);

% z-direction PML (no PML on the bottom of the grid)
psi_Exz_2 = zeros(imax-1,jmax,PML+1);
psi_Hyz_2 = zeros(imax-1,jmax,PML);

psi_Eyz_2 = zeros(imax,jmax-1,PML+1);
psi_Hxz_2 = zeros(imax,jmax-1,PML);

%**************************************************************************
% TFSF Plane Wave Source Condition
%**************************************************************************

% grid indices on which to start and end the TFSF (number of cells from the
% edge of the grid)
% the TFSF plane wave starts at k = 1 and ends at kmax-kTFSF+1
iTFSF = PML+5;
jTFSF = PML+5;
kTFSF = PML+5;  

Ez_inc = zeros(1,imax);
Hy_inc = zeros(1,imax-1); 

Ca_inc = 1.0;
Cb_inc = dt / (epso*delta);

% x-direction incident grid PML
psi_Ezx_2_inc = zeros(PML+1);
psi_Hyx_2_inc = zeros(PML);

%**************************************************************************
% Initialize the Field Components and Coefficients
%**************************************************************************

% arrays to hold the Ez and Hy fields
Ex = zeros(imax-1,jmax,kmax);
Ey = zeros(imax,jmax-1,kmax);
Ez = zeros(imax,jmax,kmax-1);

% start and end the grid on an Ez, so one less Hy value compared to Ez
Hx = zeros(imax,jmax-1,kmax-1);
Hy = zeros(imax-1,jmax,kmax-1);
Hz = zeros(imax-1,jmax-1,kmax);

% Pre-calculate the updating coefficients (so we do not have to recalculate 
% every time step). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % NEW Code segment to account for capacitance

Ca = 1.0;
Cb = dt / (epso*delta);
C_steel = 10*1E-6; % capacitance of 10uF for the lumped element at the base of crane
Cb_cap = (dt / (epso))/(1+(C_steel/(epso*delta))); % constant term for the update of Ez with capacitor at base
%steel crane
sigma = 5*1E7; % conductivity of crane
Ca_z = zeros(imax,jmax,kmax-1);
Cb_z = zeros(imax,jmax,kmax-1);
Ca_z(1:imax,1:jmax,1:kmax-1) = Ca;
Cb_z(1:imax,1:jmax,1:kmax-1) = Cb;
Ca_z(imax/2,jmax/2,2:61) = (1-((dt*sigma)/(2*epso)))/(1+((dt*sigma)/(2*epso))); % crane conuctivity update for Ca and Cb
Cb_z(imax/2,jmax/2,2:61) = (dt/(epso*delta))/(1+((sigma*dt)/(2*epso))); % crane conuctivity update for Ca and Cb
Cb_z(imax/2,jmax/2,1) = Cb_cap; % uncomment to induce the effect of the capacitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Da = 1.0;
Db = dt/(muo*delta);
current = zeros(nmax,1);

%**************************************************************************
% Start the time-stepping loop
%**************************************************************************

for n=1:nmax

%********************************
% update the magnetic field in the incident grid (Hy_inc)  
%********************************

    Hy_inc(1:imax-1) = Da*Hy_inc(1:imax-1) + ...
        Db*(Ez_inc(2:imax) - Ez_inc(1:imax-1));
    
% incident grid PML on right side (Hy components)
    for i = imax-PML:imax-1
        psi_Hyx_2_inc(imax-i) = bH(imax-i)*psi_Hyx_2_inc(imax-i) + ...
            cH(imax-i) *(Ez_inc(i+1) - Ez_inc(i))/delta;
        Hy_inc(i) = Hy_inc(i) + dt/muo*psi_Hyx_2_inc(imax-i);
    end

%********************************
% Update the Hx fields 
%********************************

    Hx(1:imax-1,1:jmax-1,1:kmax-1) = Da*Hx(1:imax-1,1:jmax-1,1:kmax-1) +...
        Db*(Ez(1:imax-1,1:jmax-1,1:kmax-1) - Ez(1:imax-1,2:jmax,1:kmax-1)...
        + Ey(1:imax-1,1:jmax-1,2:kmax) - Ey(1:imax-1,1:jmax-1,1:kmax-1));

% front TFSF for Hx
    j = jTFSF-1;
    for k = 1:kmax-kTFSF
        Hx(iTFSF:imax-iTFSF+1,j,k) = Hx(iTFSF:imax-iTFSF+1,j,k) + ...
            Db*(Ez_inc(iTFSF:imax-iTFSF+1)');
    end
   
% back TFSF for Hx
    j = jmax-jTFSF+1;
    for k = 1:kmax-kTFSF
        Hx(iTFSF:imax-iTFSF+1,j,k) = Hx(iTFSF:imax-iTFSF+1,j,k) + ...
            Db*(-Ez_inc(iTFSF:imax-iTFSF+1)');
    end

% front j-direction PML for Hx
    A = reshape(permute(psi_Hxy_1(1:imax-1,1:PML,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    B = reshape(permute(Ez(1:imax-1,1:PML,1:kmax-1)-Ez(1:imax-1,2:PML+1,...
        1:kmax-1),[1 3 2]),[(imax-1)*(kmax-1) PML]);
    C = bH(1:PML).*A' + cH(1:PML).*B'/delta;
    psi_Hxy_1(1:imax-1,1:PML,1:kmax-1) = permute(reshape(C',[(imax-1) ...
        (kmax-1) PML]),[1 3 2]);
    Hx(1:imax-1,1:PML,1:kmax-1) = Hx(1:imax-1,1:PML,1:kmax-1) + ...
        dt/muo*psi_Hxy_1(1:imax-1,1:PML,1:kmax-1);          
        
 % back j-direction PML for Hx
    A = reshape(permute(psi_Hxy_2(1:imax-1,PML:-1:1,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    B = reshape(permute(Ez(1:imax-1,jmax-PML:jmax-1,1:kmax-1)-...
        Ez(1:imax-1,jmax-PML+1:jmax,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hxy_2(1:imax-1,PML:-1:1,1:kmax-1) = permute(reshape(C',...
        [(imax-1) (kmax-1) PML]),[1 3 2]);
    Hx(1:imax-1,jmax-PML:jmax-1,1:kmax-1) = Hx(1:imax-1,...
        jmax-PML:jmax-1,1:kmax-1) + dt/muo*psi_Hxy_2(1:imax-1,PML:-1:1,...
        1:kmax-1);   
    
% top k-direction PML for Hx
    A = reshape(psi_Hxz_2(1:imax-1,1:jmax-1,PML:-1:1),...
        [(imax-1)*(jmax-1) PML]);
    B = reshape(Ey(1:imax-1,1:jmax-1,kmax-PML+1:kmax)-Ey(1:imax-1,...
        1:jmax-1,kmax-PML:kmax-1),[(imax-1)*(jmax-1) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hxz_2(1:imax-1,1:jmax-1,PML:-1:1) = reshape(C',[(imax-1) ...
        (jmax-1) PML]);
    Hx(1:imax-1,1:jmax-1,kmax-PML:kmax-1) = Hx(1:imax-1,1:jmax-1,...
        kmax-PML:kmax-1) + dt/muo*psi_Hxz_2(1:imax-1,1:jmax-1,PML:-1:1); 
           
%********************************    
% Update the Hy fields 
%********************************

    Hy(1:imax-1,2:jmax-1,1:kmax-1) = Da*Hy(1:imax-1,2:jmax-1,1:kmax-1) ...
        + Db*(Ez(2:imax,2:jmax-1,1:kmax-1) - Ez(1:imax-1,2:jmax-1,...
        1:kmax-1) + Ex(1:imax-1,2:jmax-1,1:kmax-1) - Ex(1:imax-1,...
        2:jmax-1,2:kmax));
    
% left TFSF for Hy
    i = iTFSF-1;
    Hy(i,jTFSF:jmax-jTFSF+1,1:kmax-kTFSF) = Hy(i,jTFSF:jmax-jTFSF+1,...
        1:kmax-kTFSF)+Db*(-Ez_inc(iTFSF));
        
% right TFSF for Hy
    i = imax-iTFSF+1;
    Hy(i,jTFSF:jmax-jTFSF+1,1:kmax-kTFSF) = Hy(i,jTFSF:jmax-jTFSF+1,...
        1:kmax-kTFSF)+Db*(Ez_inc(imax-iTFSF+1));
     
% left i-direction PML for Hy
    A = reshape(permute(psi_Hyx_1(1:PML,2:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-2)*(kmax-1) PML]);
    B = reshape(permute(Ez(2:PML+1,2:jmax-1,1:kmax-1)-Ez(1:PML,...
        2:jmax-1,1:kmax-1),[2 3 1]),[(jmax-2)*(kmax-1) PML]);
    C = bH(1:PML).*A' + cH(1:PML).*B'/delta;
    psi_Hyx_1(1:PML,2:jmax-1,1:kmax-1) = permute(reshape(C',[(jmax-2)...
        (kmax-1) PML]), [3 1 2]);
    Hy(1:PML,2:jmax-1,1:kmax-1) = Hy(1:PML,2:jmax-1,1:kmax-1) + ...
        dt/muo*psi_Hyx_1(1:PML,2:jmax-1,1:kmax-1);

% right i-direction PML for Hy       
    A = reshape(permute(psi_Hyx_2(PML:-1:1,2:jmax-1,1:kmax-1),...
        [2 3 1]),[(jmax-2)*(kmax-1) PML]);
    B = reshape(permute(Ez(imax-PML+1:imax,2:jmax-1,1:kmax-1)-...
        Ez(imax-PML:imax-1,2:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-2)*(kmax-1) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hyx_2(PML:-1:1,2:jmax-1,1:kmax-1) = permute(reshape(C',[(jmax-2)...
        (kmax-1) PML]), [3 1 2]);
    Hy(imax-PML:imax-1,2:jmax-1,1:kmax-1) = Hy(imax-PML:imax-1,...
        2:jmax-1,1:kmax-1) + dt/muo*psi_Hyx_2(PML:-1:1,2:jmax-1,1:kmax-1);

% top k-direction PML for Hy
    A = reshape(psi_Hyz_2(1:imax-1,2:jmax-1,PML:-1:1),...
        [(imax-1)*(jmax-2) PML]);
    B = reshape(Ex(1:imax-1,2:jmax-1,kmax-PML:kmax-1)-...
        Ex(1:imax-1,2:jmax-1,kmax-PML+1:kmax),[(imax-1)*(jmax-2) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hyz_2(1:imax-1,2:jmax-1,PML:-1:1) = reshape(C',[(imax-1) ...
        (jmax-2) PML]);
    Hy(1:imax-1,2:jmax-1,kmax-PML:kmax-1) = Hy(1:imax-1,2:jmax-1,...
        kmax-PML:kmax-1) + dt/muo*psi_Hyz_2(1:imax-1,2:jmax-1,PML:-1:1);
            
%********************************            
% Update the Hz fields 
%********************************

    Hz(1:imax-1,1:jmax-1,1:kmax-1) = Da*Hz(1:imax-1,1:jmax-1,1:kmax-1) ...
        + Db*(Ex(1:imax-1,2:jmax,1:kmax-1) - Ex(1:imax-1,1:jmax-1,...
        1:kmax-1) + Ey(1:imax-1,1:jmax-1,1:kmax-1) - Ey(2:imax,1:jmax-1,...
        1:kmax-1));
    
% left i-direction PML for Hz
    A = reshape(permute(psi_Hzx_1(1:PML,1:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-1)*(kmax-1) PML]);
    B = reshape(permute(Ey(1:PML,1:jmax-1,1:kmax-1)-Ey(2:PML+1,1:jmax-1,...
        1:kmax-1),[2 3 1]),[(jmax-1)*(kmax-1) PML]);
    C = bH(1:PML).*A' + cH(1:PML).*B'/delta;
    psi_Hzx_1(1:PML,1:jmax-1,1:kmax-1) = permute(reshape(C',[(jmax-1)...
        (kmax-1) PML]), [3 1 2]);
    Hz(1:PML,1:jmax-1,1:kmax-1) = Hz(1:PML,1:jmax-1,1:kmax-1) + ...
        dt/muo*psi_Hzx_1(1:PML,1:jmax-1,1:kmax-1);
    
% right i-direction PML for Hz 
    A = reshape(permute(psi_Hzx_2(PML:-1:1,1:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-1)*(kmax-1) PML]);
    B = reshape(permute(Ey(imax-PML:imax-1,1:jmax-1,1:kmax-1)-...
        Ey(imax-PML+1:imax,1:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-1)*(kmax-1) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hzx_2(PML:-1:1,1:jmax-1,1:kmax-1) = permute(reshape(C',[(jmax-1)...
        (kmax-1) PML]), [3 1 2]);
    Hz(imax-PML:imax-1,1:jmax-1,1:kmax-1) = Hz(imax-PML:imax-1,...
        1:jmax-1,1:kmax-1) + dt/muo*psi_Hzx_2(PML:-1:1,1:jmax-1,1:kmax-1);

% front j-direction PML for Hz
    A = reshape(permute(psi_Hzy_1(1:imax-1,1:PML,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    B = reshape(permute(Ex(1:imax-1,2:PML+1,1:kmax-1)-Ex(1:imax-1,...
        1:PML,1:kmax-1),[1 3 2]),[(imax-1)*(kmax-1) PML]);
    C = bH(1:PML).*A' + cH(1:PML).*B'/delta;
    psi_Hzy_1(1:imax-1,1:PML,1:kmax-1) = permute(reshape(C',[(imax-1)...
        (kmax-1) PML]),[1 3 2]);
    Hz(1:imax-1,1:PML,1:kmax-1) = Hz(1:imax-1,1:PML,1:kmax-1) + ...
        dt/muo*psi_Hzy_1(1:imax-1,1:PML,1:kmax-1);       
                      
% back j-direction PML for Hz
    A = reshape(permute(psi_Hzy_2(1:imax-1,PML:-1:1,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    B = reshape(permute(Ex(1:imax-1,jmax-PML+1:jmax,1:kmax-1)-...
        Ex(1:imax-1,jmax-PML:jmax-1,1:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-1) PML]);
    C = bH(PML:-1:1).*A' + cH(PML:-1:1).*B'/delta;
    psi_Hzy_2(1:imax-1,PML:-1:1,1:kmax-1) = permute(reshape(C',...
        [(imax-1) (kmax-1) PML]),[1 3 2]);
    Hz(1:imax-1,jmax-PML:jmax-1,1:kmax-1) = Hz(1:imax-1,...
        jmax-PML:jmax-1,1:kmax-1) + dt/muo*psi_Hzy_2(1:imax-1,...
        PML:-1:1,1:kmax-1);         

%********************************
% update the electric field in the incident grid (Ez_inc)  
%********************************

    Ez_inc(2:imax-1)= Ca_inc * Ez_inc(2:imax-1) + Cb_inc *...
        (Hy_inc(2:imax-1) - Hy_inc(1:imax-2));
    
    % source in the incident grid
    Ez_inc(1) = 6.11*sin(2*pi*n*freq*dt);
    
    % incident grid PML on right side (on the Ez's)
    for i = imax-PML:imax-1
        psi_Ezx_2_inc(imax-i+1) = bE(imax-i+1)*psi_Ezx_2_inc(imax-i+1)...
            + cE(imax-i+1) *(Hy_inc(i) - Hy_inc(i-1))/delta;
        Ez_inc(i) = Ez_inc(i) + dt/epso*psi_Ezx_2_inc(imax-i+1);
    end
    
%********************************            
% Update the Ex fields 
%********************************

    Ex(1:imax-1,2:jmax-1,2:kmax-1)= Ca*Ex(1:imax-1,2:jmax-1,2:kmax-1) ...
        +Cb*(Hy(1:imax-1,2:jmax-1,1:kmax-2) - Hy(1:imax-1,2:jmax-1,...
        2:kmax-1) + Hz(1:imax-1,2:jmax-1,2:kmax-1) - Hz(1:imax-1,...
        1:jmax-2,2:kmax-1));
        
% TFSF corrections on Ex
    k = kmax-kTFSF+1;
    for j = jTFSF:jmax-jTFSF+1
        Ex(iTFSF:imax-iTFSF,j,k) = Ex(iTFSF:imax-iTFSF,j,k)+...
            Cb*(-Hy_inc(iTFSF:imax-iTFSF)');
    end
    
% front j-direction PML for Ex
    A = reshape(permute(psi_Exy_1(1:imax-1,2:PML+1,2:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-2) PML]);
    B = reshape(permute(Hz(1:imax-1,2:PML+1,2:kmax-1)-Hz(1:imax-1,1:PML,...
        2:kmax-1),[1 3 2]),[(imax-1)*(kmax-2) PML]);
    C = bE(2:PML+1).*A' + cE(2:PML+1).*B'/delta;
    psi_Exy_1(1:imax-1,2:PML+1,2:kmax-1) = permute(reshape(C',[(imax-1)...
        (kmax-2) PML]),[1 3 2]);
    Ex(1:imax-1,2:PML+1,2:kmax-1) = Ex(1:imax-1,2:PML+1,2:kmax-1) + ...
        dt/epso*psi_Exy_1(1:imax-1,2:PML+1,2:kmax-1);
           
% back j-direction PML for Ex
    A = reshape(permute(psi_Exy_2(1:imax-1,PML+1:-1:2,2:kmax-1),...
        [1 3 2]),[(imax-1)*(kmax-2) PML]);
    B = reshape(permute(Hz(1:imax-1,jmax-PML:jmax-1,2:kmax-1)-...
        Hz(1:imax-1,jmax-PML-1:jmax-2,2:kmax-1),[1 3 2]),...
        [(imax-1)*(kmax-2) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Exy_2(1:imax-1,PML+1:-1:2,2:kmax-1) = permute(reshape(C',...
        [(imax-1) (kmax-2) PML]), [1 3 2]);
    Ex(1:imax-1,jmax-PML:jmax-1,2:kmax-1) = Ex(1:imax-1,jmax-PML:jmax-1,...
        2:kmax-1) + dt/epso*psi_Exy_2(1:imax-1,PML+1:-1:2,2:kmax-1);     

% top k-direction PML for Ex  
    A = reshape(psi_Exz_2(1:imax-1,2:jmax-1,PML+1:-1:2),...
        [(imax-1)*(jmax-2) PML]);
    B = reshape(Hy(1:imax-1,2:jmax-1,kmax-PML-1:kmax-2)-Hy(1:imax-1,...
        2:jmax-1,kmax-PML:kmax-1),[(imax-1)*(jmax-2) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Exz_2(1:imax-1,2:jmax-1,PML+1:-1:2) = reshape(C',[(imax-1) ...
        (jmax-2) PML]);
    Ex(1:imax-1,2:jmax-1,kmax-PML:kmax-1) = Ex(1:imax-1,2:jmax-1,...
        kmax-PML:kmax-1) + dt/epso*psi_Exz_2(1:imax-1,2:jmax-1,PML+1:-1:2);
 
%********************************            
% Update the Ey fields 
%********************************

    Ey(2:imax-1,1:jmax-1,2:kmax-1)= Ca*Ey(2:imax-1,1:jmax-1,2:kmax-1) ...
        + Cb*(Hx(2:imax-1,1:jmax-1,2:kmax-1) - Hx(2:imax-1,1:jmax-1,...
        1:kmax-2) + Hz(1:imax-2,1:jmax-1,2:kmax-1) - Hz(2:imax-1,...
        1:jmax-1,2:kmax-1));
         
% left i-direction PML for Ey   
    A = reshape(permute(psi_Eyx_1(2:PML+1,1:jmax-1,2:kmax-1),[2 3 1]),...
        [(jmax-1)*(kmax-2) PML]);
    B = reshape(permute(Hz(1:PML,1:jmax-1,2:kmax-1)-Hz(2:PML+1,1:jmax-1,...
        2:kmax-1),[2 3 1]),[(jmax-1)*(kmax-2) PML]);
    C = bE(2:PML+1).*A' + cE(2:PML+1).*B'/delta;
    psi_Eyx_1(2:PML+1,1:jmax-1,2:kmax-1) = permute(reshape(C',[(jmax-1)...
        (kmax-2) PML]),[3 1 2]);  
    Ey(2:PML+1,1:jmax-1,2:kmax-1) = Ey(2:PML+1,1:jmax-1,2:kmax-1) + ...
        dt/epso*psi_Eyx_1(2:PML+1,1:jmax-1,2:kmax-1); 
                 
% right i-direction PML for Ey
    A = reshape(permute(psi_Eyx_2(PML+1:-1:2,1:jmax-1,2:kmax-1),...
        [2 3 1]),[(jmax-1)*(kmax-2) PML]);
    B = reshape(permute(Hz(imax-PML-1:imax-2,1:jmax-1,2:kmax-1)-...
        Hz(imax-PML:imax-1,1:jmax-1,2:kmax-1),[2 3 1]),...
        [(jmax-1)*(kmax-2) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Eyx_2(PML+1:-1:2,1:jmax-1,2:kmax-1) = permute(reshape(C',...
        [(jmax-1) (kmax-2) PML]), [3 1 2]);   
    Ey(imax-PML:imax-1,1:jmax-1,2:kmax-1) = Ey(imax-PML:imax-1,1:jmax-1,...
        2:kmax-1) + dt/epso*psi_Eyx_2(PML+1:-1:2,1:jmax-1,2:kmax-1);  

% top k-direction PML for Ey
    A = reshape(psi_Eyz_2(2:imax-1,1:jmax-1,PML+1:-1:2),...
        [(imax-2)*(jmax-1) PML]);
    B = reshape(Hx(2:imax-1,1:jmax-1,kmax-PML:kmax-1)-Hx(2:imax-1,...
        1:jmax-1,kmax-PML-1:kmax-2),[(imax-2)*(jmax-1) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Eyz_2(2:imax-1,1:jmax-1,PML+1:-1:2) = reshape(C',...
        [(imax-2) (jmax-1) PML]);
    Ey(2:imax-1,1:jmax-1,kmax-PML:kmax-1) = Ey(2:imax-1,1:jmax-1,...
        kmax-PML:kmax-1) + dt/epso*psi_Eyz_2(2:imax-1,1:jmax-1,PML+1:-1:2);                
                 
%********************************            
% Update the Ez fields 
%********************************   

    Ez(2:imax-1,2:jmax-1,1:kmax-1)= Ca_z(2:imax-1,2:jmax-1,1:kmax-1).*Ez(2:imax-1,2:jmax-1,1:kmax-1) ...
        + Cb_z(2:imax-1,2:jmax-1,1:kmax-1).*(Hy(2:imax-1,2:jmax-1,...
        1:kmax-1) - Hy(1:imax-2,2:jmax-1,1:kmax-1) + Hx(2:imax-1,...
        1:jmax-2,1:kmax-1) - Hx(2:imax-1,2:jmax-1,1:kmax-1));

% TFSF corrections for Ez (left side)
    i = iTFSF;
    Ez(i,jTFSF:jmax-jTFSF+1,1:kmax-kTFSF) = Ez(i,jTFSF:jmax-jTFSF+1,...
        1:kmax-kTFSF)+Cb*(-Hy_inc(i-1));

% TFSF corrections for Ez (right side)
    i = imax-iTFSF+1;
    Ez(i,jTFSF:jmax-jTFSF+1,1:kmax-kTFSF) = Ez(i,jTFSF:jmax-jTFSF+1,...
        1:kmax-kTFSF)+Cb*(Hy_inc(i));
        
% left i-direction PML on Ez
    A = reshape(permute(psi_Ezx_1(2:PML+1,2:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-2)*(kmax-1) PML]);
    B = reshape(permute(Hy(2:PML+1,2:jmax-1,1:kmax-1)-Hy(1:PML,2:jmax-1,...
        1:kmax-1),[2 3 1]),[(jmax-2)*(kmax-1) PML]);
    C = bE(2:PML+1).*A' + cE(2:PML+1).*B'/delta;
    psi_Ezx_1(2:PML+1,2:jmax-1,1:kmax-1) = permute(reshape(C',[(jmax-2)...
        (kmax-1) PML]), [3 1 2]);
    Ez(2:PML+1,2:jmax-1,1:kmax-1) = Ez(2:PML+1,2:jmax-1,1:kmax-1) + ...
        dt/epso*psi_Ezx_1(2:PML+1,2:jmax-1,1:kmax-1);

% right i-direction PML on Ez
    A = reshape(permute(psi_Ezx_2(PML+1:-1:2,2:jmax-1,1:kmax-1),...
        [2 3 1]),[(jmax-2)*(kmax-1) PML]);
    B = reshape(permute(Hy(imax-PML:imax-1,2:jmax-1,1:kmax-1)-...
        Hy(imax-PML-1:imax-2,2:jmax-1,1:kmax-1),[2 3 1]),...
        [(jmax-2)*(kmax-1) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Ezx_2(PML+1:-1:2,2:jmax-1,1:kmax-1) = permute(reshape(C',...
        [(jmax-2) (kmax-1) PML]), [3 1 2]);
    Ez(imax-PML:imax-1,2:jmax-1,1:kmax-1) = Ez(imax-PML:imax-1,2:jmax-1,...
        1:kmax-1) + dt/epso*psi_Ezx_2(PML+1:-1:2,2:jmax-1,1:kmax-1);
    
% front j-direction PML on Ez
    A = reshape(permute(psi_Ezy_1(2:imax-1,2:PML+1,1:kmax-1),[1 3 2]),...
        [(imax-2)*(kmax-1) PML]);
    B = reshape(permute(Hx(2:imax-1,1:PML,1:kmax-1)-Hx(2:imax-1,2:PML+1,...
        1:kmax-1),[1 3 2]),[(imax-2)*(kmax-1) PML]);
    C = bE(2:PML+1).*A' + cE(2:PML+1).*B'/delta;
    psi_Ezy_1(2:imax-1,2:PML+1,1:kmax-1) = permute(reshape(C',[(imax-2)...
        (kmax-1) PML]),[1 3 2]);
    Ez(2:imax-1,2:PML+1,1:kmax-1) = Ez(2:imax-1,2:PML+1,1:kmax-1) + ...
        dt/epso*psi_Ezy_1(2:imax-1,2:PML+1,1:kmax-1);
    
% back j-direction PML on Ez
    A = reshape(permute(psi_Ezy_2(2:imax-1,PML+1:-1:2,1:kmax-1),...
        [1 3 2]),[(imax-2)*(kmax-1) PML]);
    B = reshape(permute(Hx(2:imax-1,jmax-PML-1:jmax-2,1:kmax-1)-...
        Hx(2:imax-1,jmax-PML:jmax-1,1:kmax-1),[1 3 2]),...
        [(imax-2)*(kmax-1) PML]);
    C = bE(PML+1:-1:2).*A' + cE(PML+1:-1:2).*B'/delta;
    psi_Ezy_2(2:imax-1,PML+1:-1:2,1:kmax-1) = permute(reshape(C',...
        [(imax-2) (kmax-1) PML]),[1 3 2]);
    Ez(2:imax-1,jmax-PML:jmax-1,1:kmax-1) = Ez(2:imax-1,jmax-PML:jmax-1,...
        1:kmax-1) + dt/epso*psi_Ezy_2(2:imax-1,PML+1:-1:2,1:kmax-1);
 
%********************************            
% Output and Plotting
%********************************   
 current(n) = delta*(Hy(imax/2,jmax/2,31)- Hy(imax/2-1,jmax/2,31) - Hx(imax/2,jmax/2,31) + Hz(imax/2,jmax/2-1,31));


 

 if (mod(n,500) == 0) 
        
    % Plot the Ez values in the x-y plane at k = kmax/2 during time stepping
    figure(1);
    pcolor((Ez(:,:,kmax/2)'));
    shading flat;
    colorbar;
    title('x-y plane of Ez fields at k = kmax/2', 'FontSize',16);
    xlabel('i index', 'FontSize',14);
    ylabel('j index', 'FontSize',14);
    set(gca,'FontSize',14);

    % Plot the Ez values in the x-z plane at j = jmax/2 during time stepping
    figure(2);
    for i = 1:imax
        for k = 1:kmax-1
            Ezplot(i,k) = Ez(i,jmax/2,k);
        end
    end
    pcolor(Ezplot(:,:)'); 
    shading flat;
    colorbar;
    title('x-z plane of Ez fields at j = jmax/2', 'FontSize',16);
    xlabel('i index (this side of the grid is the ground)', 'FontSize',14);
    ylabel('k index', 'FontSize',14);
    set(gca,'FontSize',14);

    pause(0.001);
    
    % plot the current at the midpoint of the crane
    
    
 end

%********************************            
% End the time-stepping loop
%********************************

end

figure(3);
plot(current,'LineWidth',2);
xlabel('Time Steps n', 'FontSize',14);
ylabel('Current(A/m)', 'FontSize',14);
title('Current Monitor, Crane Only', 'FontSize',14);


