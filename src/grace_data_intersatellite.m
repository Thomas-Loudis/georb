function [intersat_struct] = grace_data_intersatellite (cfg_fname, intersat_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: grace_data_intersatellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read intersatellite ranging data from GRACE mission 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  27 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 03/05/2025  Thomas Loudis Papanikolaou
%             Code minor modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intersatellite ranging data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'intersat_ranging_data_simul';
[param_value] = read_param_cfg(cfg_fname,param_keyword);
intersat_ranging_data_simul = param_value;
test_intersat_ranging_data_simul = strcmp(intersat_ranging_data_simul,'y');

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

% KBR or LRI data file name (GNV format) or Simulated Ranging Data file name (georb format)
param_keyword = intersat_type;
[param_value] = read_param_cfg(cfg_fname,param_keyword);
intersat_ranging_data_filename = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated Intersatellite Ranging data (in georb format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test_intersat_ranging_data_simul == 1
    % GEORB Intersatellite ranging data
    [intersat_data_matrix, header_data_matrix] = read_georb_data(intersat_ranging_data_filename);

    [sz1, sz2] = size(intersat_data_matrix);    
    intersat_ranging = zeros(sz1,4);

    % MJD in TT including fraction of the day
    intersat_ranging(:,1) = intersat_data_matrix(:,1) + intersat_data_matrix(:,2) / 86400;
    % Range
    intersat_ranging(:,2) = intersat_data_matrix(:,3);
    % Range-Rate
    intersat_ranging(:,3) = intersat_data_matrix(:,4);
    % Range-acceleration 
    intersat_ranging(:,4) = intersat_data_matrix(:,5);

    biasrange = [intersat_ranging(:,1) intersat_ranging(:,2)];
    rangerate = [intersat_ranging(:,1) intersat_ranging(:,3)];
    rangeaccl = [intersat_ranging(:,1) intersat_ranging(:,4)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE missions KBR/LRI reading and preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kbr_data_fname = intersat_ranging_data_filename;

% Read intersatellite ranging data (corrections are added during reading via function grace_kbr1b.m)
[kbr1b,biasrange,rangerate,rangeaccl] = grace_kbr1b(kbr_data_fname);

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
end    


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
