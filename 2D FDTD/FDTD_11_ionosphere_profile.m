sigma = zeros(kmax-1,1);
electron_charge = -1.6E-19;  % electron charge
electron_mass = 9.1E-31;     % electron mass
sharpness = 0.5;             % sharpness of the ionosphere profile
reflect_altitude = 75.0;     % reflection height of the ionosphere

 for k = 1:kmax-1
     collision_freq = 1.816e11 * exp(-0.15 * k*delta/1000);  
     electron_density = 1.43E13 * exp(-0.15 * reflect_altitude) * ...
         exp((sharpness-0.15)*(k*delta/1000 - reflect_altitude));
     sigma(k) = (electron_density * electron_charge^2) / ...
         (electron_mass * collision_freq);
 end
 