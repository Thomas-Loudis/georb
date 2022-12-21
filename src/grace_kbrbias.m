function [KBRbias,nonbiasrange,resrange,rms_resrange,resrangerate,rms_resrangerate] = grace_kbrbias(orbA,orbB,biasrange,rangerate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_kbrbias:  GRACE K-band range (KBR) bias estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  KBR bias estimation based on orbit data
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit differences: XYZ
[dorb,rms_orb] = compstat(orbA,orbB);
% Satellites distance via orbit differences XYZ
[sz1 sz2] = size(dorb);
dr = zeros(sz1,2);
dv = zeros(sz1,2);
for i = 1 : sz1
    dr(i,:) = [dorb(i,1) sqrt(dorb(i,2)^2 + dorb(i,3)^2 + dorb(i,4)^2)];
    dv(i,:) = [dorb(i,1) sqrt(dorb(i,5)^2 + dorb(i,6)^2 + dorb(i,7)^2)];
end
clear sz1 sz2 i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR bias estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find common KBR1B biased ranges and orbits distance
[drkbr,rms_drkbr,rkbr] = compstat(dr,biasrange);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biased range data refer to range plus an unknown bias
% KBRbias = 1/n * Sum( biasedrange - dr(XYZ) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(drkbr);
sum_drange = 0;
for i = 1 : sz1
    sum_drange = sum_drange + drkbr(i,2);
end
KBRbias = (1/sz1) * (sum_drange);
clear i sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-biased ranges: Remove bias from biased range
[sz1 sz2] = size(biasrange);
nonbiasrange = zeros(sz1,2);
for i = 1 : sz1
    nonbiasrange(i,:) = [biasrange(i,1) (biasrange(i,2)-KBRbias)];
end
clear i sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting (residuals) of used orbits to nonbiasrange and rangerate
[resrange,rms_resrange] = compstat(dr,nonbiasrange);
[resrangerate,rms_resrangerate] = compstat(dv,rangerate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


