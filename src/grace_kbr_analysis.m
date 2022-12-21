function [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, infilenameGA,infilenameGB,orbcA,orbcB,intersat_obs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%         KBR analysis / residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - infilename  : *.in input file name
% - infilename2 : modified *.in input filename that changes are written 
% - infilename3 : name of the file with the modified prmlines to be written 
%
% Output arguments:
% - orbc:         State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:            Epoch in seconds of TT (Terrestrial Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         July 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 08/07/2022   Thomas Loudis Papanikolaou 
%              Code upgrade based on the new orbit configuration format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR/LRI reading and preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case: KBR
test_intersat_obs_type = strcmp(intersat_obs,'intersat_KBR');
if test_intersat_obs_type == 1
param_keyword = 'KBR_data';
end
% Case: LRI
test_intersat_obs_type = strcmp(intersat_obs,'intersat_LRI');
if test_intersat_obs_type == 1
param_keyword = 'LRI_data';
end

% KBR or LRI data file name
[param_value] = read_param_cfg(infilenameGA,param_keyword);
kbr_data_fname = param_value;

% Read intersatellite ranging data (corrections are added during reading via function grace_kbr1b.m)
tstop = 0;
[kbr1b,biasrange,rangerate,rangeaccl] = grace_kbr1b(kbr_data_fname, tstop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change format : remove parameter t(TT) sec since 0h (2nd column)
kbr1b = [kbr1b(:,1) kbr1b(:,3:end)];
biasrange = [biasrange(:,1) biasrange(:,3:end)]; 
rangerate = [rangerate(:,1) rangerate(:,3:end)]; 
rangeaccl = [rangeaccl(:,1) rangeaccl(:,3:end)]; 
% Time conversion from GPS time to Terrestrial Time scale
[sz1 sz2] = size(kbr1b);
for i = 1 : sz1
    mjdgps = kbr1b(i,1);
    [tgps,D,M,Y] = MJD_inv(mjdgps);
    [tutc,tTT] = time_scales_GPS(tgps,mjdgps);
    mjdTT = mjdgps + (tTT-tgps)/60/60/24;
    kbr1b(i,1) = mjdTT;
    biasrange(i,1) = mjdTT;
    rangerate(i,1) = mjdTT;
    rangeaccl(i,1) = mjdTT;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR bias estimation (based on GNV1B orbits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KBRbias,nonbiasrange,resrangeGNV1B,rms_resrangeGNV1B,resrangerateGNV1B,rms_resrangerateGNV1B] = grace_kbrbias(orbcA,orbcB,biasrange,rangerate);
KBRbias_orbc  = KBRbias;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR analysis / residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[resrange,rms_resrange,resrangerate,rms_resrangerate,dresrange,dresrangerate,rms_dresrange,rms_dresrangerate] = grace_kbrdres(orbcA,orbcB,nonbiasrange,rangerate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
