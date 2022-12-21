function [resrange,rms_resrange,resrangerate,rms_resrangerate,dresrange,dresrangerate,rms_dresrange,rms_dresrangerate] = grace_kbrdres(orbA,orbB,range,rangerate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_kbrdres:  GRACE KBR residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  KBR residuals from estimated GRACE orbits. Range and range-rate
%  residuals are computed as well as delta_range and delta_range-rate based
%  on epoch differences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:   KBR1b data file's name
%
% Output arguments:
% - kbr1b:      KBR data array
%   kbr1b = [MJDgps tgps KBR1B(2:end)]
% - biasrange:  Corrected biased range data
%   biasrange = [MJDgps tgps biased_range]
% - rangerate:  Corrected range rate data
%   rangerate = [MJDgps tgps range_rate]
% - rangeaccl:  Corrected range acceleration data
%   rangeaccl = [MJDgps tgps range_acceleration]
%
%   MJDgps:     MJD in GPS time including fraction of the day
%   tgps:       Seconds since 0h in GPS time scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 19/07/2012   computation of range-rate residuals has been revised
%              dv array has been modified 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit differences: XYZ
[dstn,rms_stn,dorb,rms_orb,delta_kepler,rms_kepler,delta_Vstn rms_Vstn] = orbital_pert(orbA,orbB,0);
%[dorb,rms_orb] = compstat(orbA,orbB);

% Satellites distance via orbit differences XYZ
[sz1 sz2] = size(dorb);
dr = zeros(sz1,2);
for i = 1 : sz1
    dr(i,:) = [dorb(i,1) sqrt(dorb(i,2)^2 + dorb(i,3)^2 + dorb(i,4)^2)];
    rab_vec = [dorb(i,2) dorb(i,3) dorb(i,4)]';
    rab_magn = sqrt(dorb(i,2)^2 + dorb(i,3)^2 + dorb(i,4)^2);
    % Line-Of-Sight vector
    eab_vec = (1 / rab_magn) * rab_vec;
    vab_vec = [dorb(i,5) dorb(i,6) dorb(i,7)]';
    % range-rate
    %rangerate_orbits = vab_vec' * eab_vec;
    dv(i,:) = [ dorb(i,1)  (vab_vec' * eab_vec) ];
end
clear sz1 sz2 i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR residuals
% Residual = KBRrange - SXYZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[resrange,rms_resrange,range_common] = compstat(dr,range);
[resrangerate,rms_resrangerate,resrangerate_common] = compstat(dv,rangerate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR residuals (epoch differences)
% bias is eliminated by epoch differences
% dResidual = residual(i) - residual(i-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(range_common);
j = 1;
for i = 1 : sz1
    if i > 1
        dresrange(j,:) =[ (range_common(i,1) - range_common(i-1,1))  (range_common(i,3) - range_common(i,2)) - (range_common(i-1,3) - range_common(i-1,2)) ];
        dresrangerate(j,:) = [ (resrangerate_common(i,1) - resrangerate_common(i-1,1))  (resrangerate_common(i,3) - resrangerate_common(i,2)) - (resrangerate_common(i-1,3) - resrangerate_common(i-1,2)) ];    
        j = j + 1;
    end
end
clear sz1 sz2 i j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
[sz1 sz2] = size(dresrange);
j = 1;
for i = 2 : sz2
    rms_dresrange(1,j) = rms(dresrange(:,i));
    rms_dresrangerate(1,j) = rms(dresrangerate(:,i));
    j = j + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
