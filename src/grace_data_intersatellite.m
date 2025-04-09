function [intersat_struct] = grace_data_intersatellite (cfg_fname, intersat_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: grace_data_intersatellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read intersatellite ranging data from GRACE mission 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  27 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE Intersatellite ranging data :: KBR & LRI data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR/LRI reading and preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Case: KBR
% test_intersat_obs_type = strcmp(intersat_obs,'intersat_KBR');
% if test_intersat_obs_type == 1
% param_keyword = 'KBR_data';
% end
% % Case: LRI
% test_intersat_obs_type = strcmp(intersat_obs,'intersat_LRI');
% if test_intersat_obs_type == 1
% param_keyword = 'LRI_data';
% end

% intersat_type = LRI_data;
% intersat_type = KBR_data;
% param_keyword = 'KBR_data';

param_keyword = intersat_type;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR or LRI data file name
[param_value] = read_param_cfg(cfg_fname,param_keyword);
kbr_data_fname = param_value;

% Read intersatellite ranging data (corrections are added during reading via function grace_kbr1b.m)
[kbr1b,biasrange,rangerate,rangeaccl] = grace_kbr1b(kbr_data_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format modification :: remove parameter t(TT) sec since 0h (2nd column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kbr1b = [kbr1b(:,1) kbr1b(:,3:end)];
biasrange = [biasrange(:,1) biasrange(:,3:end)]; 
rangerate = [rangerate(:,1) rangerate(:,3:end)]; 
rangeaccl = [rangeaccl(:,1) rangeaccl(:,3:end)]; 

% Time conversion from GPS time to Terrestrial Time scale
[sz1 sz2] = size(kbr1b);
for i = 1 : sz1
    mjdgps = kbr1b(i,1);
    [tgps,D,M,Y] = MJD_inv(mjdgps);
    % [tutc,tTT] = time_scales_GPS(tgps,mjdgps);
    tTT = tgps + 51.184;
    mjdTT = mjdgps + (tTT-tgps)/60/60/24;
    kbr1b(i,1) = mjdTT;
    biasrange(i,1) = mjdTT;
    rangerate(i,1) = mjdTT;
    rangeaccl(i,1) = mjdTT;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE intersatellite data use in parameter estimation y/n 
% intersat_struct.param_estim_yn = acc_data_yn;
% Intersatellite Range data 
intersat_struct.range = biasrange;
% Intersatellite Range-Rate data 
intersat_struct.rangerate = rangerate;
% Intersatellite Range-Acceleration data 
intersat_struct.rangeacceleration = rangeaccl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
