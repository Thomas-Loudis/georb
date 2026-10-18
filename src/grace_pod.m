function [grace_pod_struct, out_dir_name] = grace_pod(grace_pod_struct, intersat_ranging_flag, intersat_obs_flag, write_data_01) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function: grace_pod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%   Orbit determination of GRACE satellites and inter-satellite ranging functionals and residuals      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - grace_pod_struct        : Overall structure array for orbit models and matrices for both satellites and intersatellite ranging
% - intersat_ranging_flag   : 'y' or 'n' InterSatellite Ranging Residuals computation based on initial observations (without empirical parameters)
% - intersat_obs_flag       : 'LRI' or 'KBR' for include/exclude LRI ranging residuals computation (GRACE-FO case)   
% - write_data_01           : 0/1 Write matrices to Data files 
% 
% Output arguments: 
% - grace_pod_struct        : Overall structure array with updated values 
% - out_dir_name            : Output Data folder name  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                              20 April 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 28/09/2026, Thomas Papanikolaou 
%             Source Code minor modification to Input arguments   
%-------------------------------------------------------------------------- 

 
%-------------------------------------------------------------------------- 
% POD of GRACE satellites 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Read Structure Arrays from the main grace_pod structure array 
%-------------------------------------------------------------------------- 
% Orbit modelling strucure arrays
orbit_model_GRACE1 = grace_pod_struct.grace1_pod.orbit_model; 
orbit_model_GRACE2 = grace_pod_struct.grace2_pod.orbit_model; 

% Orbit configuration strucutre arrays 
orbit_config_struct_GRACE1 = orbit_model_GRACE1.orbit_config; 
orbit_config_struct_GRACE2 = orbit_model_GRACE2.orbit_config; 

% Orbit matrices structure arrays
orbit_matrices_GRACE1 = grace_pod_struct.grace1_pod.orbit_matrices; 
orbit_matrices_GRACE2 = grace_pod_struct.grace2_pod.orbit_matrices; 
%-------------------------------------------------------------------------- 

write_data = 0; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit Determination of GRACE-A/-C 
[orbit_config_G1, orbit_model_GRACE1, orbit_matrices_GRACE1, OUT_data_foldername_GRACE1] = orbit_object(orbit_config_struct_GRACE1, write_data, orbit_model_GRACE1, orbit_matrices_GRACE1);  
% Update main structure array to Output
grace_pod_struct.grace1_pod.orbit_config    = orbit_config_G1; 
grace_pod_struct.grace1_pod.orbit_model     = orbit_model_GRACE1; 
grace_pod_struct.grace1_pod.orbit_matrices  = orbit_matrices_GRACE1; 
grace_pod_struct.grace1_pod.data_output_foldername = OUT_data_foldername_GRACE1; 
% grace_pod_struct.grace1_pod   = grace1_pod; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit Determination of GRACE-B/-D 
[orbit_config_G2, orbit_model_GRACE2, orbit_matrices_GRACE2, OUT_data_foldername_GRACE2] = orbit_object(orbit_config_struct_GRACE2, write_data, orbit_model_GRACE2, orbit_matrices_GRACE2); 
% Update main structure array to Output
grace_pod_struct.grace2_pod.orbit_config    = orbit_config_G2; 
grace_pod_struct.grace2_pod.orbit_model     = orbit_model_GRACE2; 
grace_pod_struct.grace2_pod.orbit_matrices  = orbit_matrices_GRACE2; 
grace_pod_struct.grace2_pod.data_output_foldername = OUT_data_foldername_GRACE2; 
% grace_pod_struct.grace2_pod   = grace2_pod; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% orbit.matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
orbcGA = orbit_matrices_GRACE1.orbit_crf; 
orbcGB = orbit_matrices_GRACE2.orbit_crf; 
 
% rms_obs_GA    = orbit_matrices_GRACE1.rms_observation_residuals; 
% rms_orbitalGA = orbit_matrices_GRACE1.rms_orbital; 
 
% rms_obs_GB    = orbit_matrices_GRACE2.rms_observation_residuals; 
% rms_orbitalGB = orbit_matrices_GRACE2.rms_orbital; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Intersatellite Ranging KBR/LRI data residuals  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inter-Satellite Ranging residuals computation test-flag
test_intersat_ranging    = strcmp(intersat_ranging_flag,'y'); 
if test_intersat_ranging == 1 
 
%-------------------------------------------------------------------------- 
% K-Band Ranging (KBR) data residuals 
intersat_obs = 'intersat_KBR'; 
[biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ... 
 resrange, resrangerate, dresrange, dresrangerate, ... 
 rms_resrange, KBR_rms_resrangerate, rms_dresrange, rms_dresrangerate, ... 
 intersat_ranging_functionals] ... 
 = grace_kbr_analysis(orbcGA, orbcGB, intersat_obs, orbit_model_GRACE1); 
KBR_rms_res_rangerate           = KBR_rms_resrangerate; 
KBR_rms_res_range               = rms_resrange; 
KBR_range_residuals             = resrange; 
KBR_rangerate_residuals         = resrangerate; 
KBR_intersat_observation_data   = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)]; 
%--------------------------------------------------------------------------  
% Laser Ranging Interferometry LRI data residuals       
LRI_obs_analysis_01 = 0; 
param_keyword = 'LRI'; 
[test_intersat_obs_flag] = strcmp(intersat_obs_flag,param_keyword); 
if test_intersat_obs_flag == 1 
    LRI_obs_analysis_01 = 1; 
end 

if LRI_obs_analysis_01 == 1 
    intersat_obs = 'intersat_LRI'; 
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ... 
     resrange, resrangerate, dresrange, dresrangerate, ... 
     rms_resrange, LRI_rms_resrangerate, rms_dresrange, rms_dresrangerate, ... 
     intersat_ranging_functionals] ... 
     = grace_kbr_analysis(orbcGA, orbcGB, intersat_obs, orbit_model_GRACE1);       
    LRI_rms_res_rangerate           = LRI_rms_resrangerate; 
    LRI_rms_res_range               = rms_resrange; 
    LRI_range_residuals             = resrange; 
    LRI_rangerate_residuals         = resrangerate;     
    LRI_intersat_observation_data   = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)];   
else 
    LRI_rms_res_rangerate = 0; 
    LRI_rms_res_range     = 0;     
end 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Inter-Satellite Observations matrix    
%-------------------------------------------------------------------------- 
if LRI_obs_analysis_01 == 1 
% LRI observations 
intersat_observation_data = LRI_intersat_observation_data; 
else 
% KBR observations 
intersat_observation_data = KBR_intersat_observation_data; 
end 
intersat_range_residuals = resrange; 
intersat_rangerate_residuals = resrangerate; 
%-------------------------------------------------------------------------- 
 
if 1>0 
fprintf('%s \n', 'Intersatellite Ranging residuals:'); 
if LRI_obs_analysis_01 == 1 
fprintf('%s%21.6f', 'LRI range residuals: RMS (mm):          ',LRI_rms_res_range * 10^3); 
fprintf('\n'); 
fprintf('%s%17.6f', 'LRI range-rate residuals: RMS (μm/sec): ',LRI_rms_res_rangerate * 10^6); 
fprintf('\n'); 
end 
fprintf('%s%17.6f', 'KBR range residuals: RMS (mm):          ',KBR_rms_res_range * 10^3); 
fprintf('\n'); 
fprintf('%s%17.6f', 'KBR range-rate residuals: RMS (μm/sec): ',KBR_rms_res_rangerate * 10^6); 
fprintf('\n'); 
fprintf('\n') 
end 
 
%-------------------------------------------------------------------------- 
% Intersatellite ranging Output structure array 
intersat_pod.intersat_observation_data    = intersat_observation_data; 
intersat_pod.intersat_range_residuals     = intersat_range_residuals; 
intersat_pod.intersat_rangerate_residuals = intersat_rangerate_residuals; 
 
intersat_pod.KBR_rms_rangerate_residuals  = KBR_rms_res_rangerate; 
intersat_pod.KBR_rms_range_residuals      = KBR_rms_res_range; 
 
intersat_pod.LRI_rms_rangerate_residuals  = LRI_rms_res_rangerate; 
intersat_pod.LRI_rms_range_residuals      = LRI_rms_res_range; 
 
% Overall Output structure array 
grace_pod_struct.intersat_pod = intersat_pod; 
grace_pod_struct.intersat_ranging_functionals = intersat_ranging_functionals; 
%-------------------------------------------------------------------------- 

else

%-------------------------------------------------------------------------- 
% Intersatellite Ranging functionals computed based on provided orbits  
[intersat_range, intersat_rangerate, intersat_ranging_functionals] = intersat_ranging(orbcGA, orbcGB); 
% Update values to grace_pod_struct 
grace_pod_struct.intersat_ranging_functionals = intersat_ranging_functionals; 
%-------------------------------------------------------------------------- 
intersat_pod.update = 'no';
 
end
%-------------------------------------------------------------------------- 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Write computed GRACE data to directory 
% Wirte data to output folder in georb format 
write_data = write_data_01; 
if write_data == 1 
neq_flag = 0; 
[OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag, neq_flag); 
out_dir_name = OUT_fname_mission_mjd; 
else 
out_dir_name = '';    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
