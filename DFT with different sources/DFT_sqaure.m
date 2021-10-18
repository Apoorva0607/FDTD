% DFT code for sqaure wave source
clc;
clear;
%% Initialization
delta = 4.3*1e-3; % dx
c = 2.99792458*1e+08;% m/s speed of light
S = 0.99;
dt = (S*delta)/c;
% fstart and fend were calculated based on the ISM band
fstart = input('Enter your start frequency'); % fstart = 820*1e6;
fend = input('Enter your end frequency');% fend = 980*1e6;
nfreq = 1000; % no. of frequency steps
df = (fend-fstart)/(nfreq-1);
omega = zeros(nfreq,1); % defining omega
for f = 1:nfreq
    omega(f) = 2*pi*fstart + 2*pi*df*(f-1);
end

ftr = zeros(nfreq,1);
fti = zeros(nfreq,1);
mag = zeros(nfreq,1);

%% source - sqaure pulse or gaussian pulse
n = input('Select 1 for square wave and 2 for gaussian wave: ');
switch n
    case 1
        nmax = 50;
        time_waveform = zeros(nmax,1);
        time_waveform(1:20) = 1;
    case 2
        nmax = 20000;
        time_waveform = zeros(nmax,1);
        f0 = 915*1e6;
        thalf = 2/(pi*13*1e6);
        t0 = 3*thalf;
        for i = 1:nmax
            time_waveform(i) = sin(2*pi*f0*(i*dt-t0))*exp(-((i*dt-t0)/(thalf))^2);
        end        
end

figure(1);
plot(time_waveform,'LineWidth',2);
title('Source waveform - sqaure pulse', 'LineWidth', 14);
xlabel('time(s)','LineWidth', 14);
ylabel('amplitude','LineWidth', 14);

%% numerical DFT of the source pulse
for n = 1:nmax
    for f = 1:nfreq
        ftr(f) = ftr(f) + time_waveform(n)*cos(omega(f)*n*dt);
        fti(f) = fti(f) + time_waveform(n)*sin(omega(f)*n*dt);
        mag(f) = sqrt(ftr(f)^2 + fti(f)^2);
    end
end


figure(2); plot(omega./(2*pi),mag./(max(mag)),'LineWidth', 2);
title('source spectrum','LineWidth', 14);
xlabel('Frequency','LineWidth', 14)
ylabel('Magnitude', 'LineWidth', 14);

