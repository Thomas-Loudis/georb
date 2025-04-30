function [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(orbcA,orbcB,intersat_obs, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Intersatellite Ranging data KBR/LRI residuals
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
% 07/04/2025   Thomas Loudis Papanikolaou
%              Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Force model matrix 
KBR_intersat_glob = orbit_model_struct.KBR_intersat_struct;
LRI_intersat_glob = orbit_model_struct.LRI_intersat_struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case: KBR
test_intersat_obs_type = strcmp(intersat_obs,'intersat_KBR');
if test_intersat_obs_type == 1
param_keyword = 'KBR_data';
intersat_struct = KBR_intersat_glob;
end

% Case: LRI
test_intersat_obs_type = strcmp(intersat_obs,'intersat_LRI');
if test_intersat_obs_type == 1
param_keyword = 'LRI_data';
intersat_struct = LRI_intersat_glob;
end

biasrange = intersat_struct.range;
rangerate = intersat_struct.rangerate;
rangeaccl = intersat_struct.rangeacceleration;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR bias estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KBRbias,nonbiasrange,resrangeGNV1B,rms_resrangeGNV1B,resrangerateGNV1B,rms_resrangerateGNV1B] = grace_kbrbias(orbcA,orbcB,biasrange,rangerate);
KBRbias_orbc  = KBRbias;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[resrange,rms_resrange,resrangerate,rms_resrangerate,dresrange,dresrangerate,rms_dresrange,rms_dresrangerate] = grace_kbrdres(orbcA,orbcB,nonbiasrange,rangerate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
