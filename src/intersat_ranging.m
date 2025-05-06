function [intersat_range, intersat_rangerate, intersat_ranging_struct] = intersat_ranging(orbA, orbB)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersat_ranging:  Intersatellite Ranging data based on orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Intersatellite Ranging functionals (range and range-rate) computed based
%  on probided input orbits.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:   KBR1b data file's name
%
% Output arguments:
% - intersat_range:         Range matrix
%   intersat_range = [MJD range]
% - intersat_rangerate:     Range rate matrix
%   intersat_rangerate = [MJD range_rate]
% - intersat_rangeaccel:    Range acceleration data
%   intersat_rangeaccel = [MJD range_acceleration]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                                 1 May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit differences: XYZ
[dstn,rms_stn,dorb,rms_orb,delta_kepler,rms_kepler,delta_Vstn rms_Vstn] = orbital_pert(orbA,orbB,0);
%[dorb,rms_orb] = compstat(orbA,orbB);

% Satellites distance via orbit differences XYZ
[sz1, sz2] = size(dorb);
dr = zeros(sz1,2); 
dv = zeros(sz1,2);
for i = 1 : sz1
    dr(i,:) = [dorb(i,1) sqrt(dorb(i,2)^2 + dorb(i,3)^2 + dorb(i,4)^2)];
    rab_vec = [dorb(i,2) dorb(i,3) dorb(i,4)]';
    rab_magn = sqrt(dorb(i,2)^2 + dorb(i,3)^2 + dorb(i,4)^2);

    % Line-Of-Sight vector
    eab_vec = (1 / rab_magn) * rab_vec;
    vab_vec = [dorb(i,5) dorb(i,6) dorb(i,7)]';
    
    % range-rate
    % rangerate_orbits = vab_vec' * eab_vec;
    dv(i,:) = [ dorb(i,1)  (vab_vec' * eab_vec) ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite ranging computed based on orbits :: intersat_ranging_orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intersat_range      = dr;
intersat_rangerate  = dv;
intersat_ranging_data = [dr(:,1) dr(:,2) dv(:,2) zeros(sz1,1)];

intersat_ranging_struct.range       = dr;
intersat_ranging_struct.rangerate   = dv;
intersat_ranging_struct.rangingdata = intersat_ranging_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

