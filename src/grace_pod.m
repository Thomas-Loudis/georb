function [grace_pod_struct, grace1_pod, grace2_pod, intersat_pod, out_dir_name] = grace_pod(orbit_model_GRACE1, orbit_model_GRACE2, intersat_obs_flag, write_data_01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: grace_pod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_model_struct  : Orbit model structure array 
% - ic_data_struct      : Initial Conditions structure array 
%
% Output arguments:
% - out_dir_name        : Output Data folder name 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                              20 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD of GRACE satellites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit configuration strucutre arrays
orbit_config_struct_GRACE1 = orbit_model_GRACE1.orbit_config;
orbit_config_struct_GRACE2 = orbit_model_GRACE2.orbit_config;

write_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination of GRACE-A/-C
[orbit_config_G1, orbit_model_GRACE1, orbit_matrices_GRACE1, OUT_data_foldername_GRACE1] = orbit_object(orbit_config_struct_GRACE1, write_data, orbit_model_GRACE1); 
grace1_pod.orbit_config = orbit_config_G1;
grace1_pod.orbit_model = orbit_model_GRACE1;
grace1_pod.orbit_matrices = orbit_matrices_GRACE1;
grace1_pod.data_output_foldername = OUT_data_foldername_GRACE1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination of GRACE-B/-D
[orbit_config_G2, orbit_model_GRACE2, orbit_matrices_GRACE2, OUT_data_foldername_GRACE2] = orbit_object(orbit_config_struct_GRACE2, write_data, orbit_model_GRACE2);
grace2_pod.orbit_config = orbit_config_G2;
grace2_pod.orbit_model = orbit_model_GRACE2;
grace2_pod.orbit_matrices = orbit_matrices_GRACE2;
grace2_pod.data_output_foldername = OUT_data_foldername_GRACE2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orbit.matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orbcGA = orbit_matrices_GRACE1.orbit_crf;
orbcGB = orbit_matrices_GRACE2.orbit_crf;

rms_obs_GA    = orbit_matrices_GRACE1.rms_observation_residuals;
rms_orbitalGA = orbit_matrices_GRACE1.rms_orbital;

rms_obs_GB    = orbit_matrices_GRACE2.rms_observation_residuals;
rms_orbitalGB = orbit_matrices_GRACE2.rms_orbital;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intersatellite Ranging KBR/LRI data residuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRI_obs_analysis_01 = 0;
param_keyword = 'LRI';
[intersat_obs_flag] = strcmp(intersat_obs_flag,param_keyword);
if intersat_obs_flag == 1
    LRI_obs_analysis_01 = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laser Ranging Interferometry LRI data residuals      
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Satellite Observations matrix   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LRI_obs_analysis_01 == 1
% LRI observations
intersat_observation_data = LRI_intersat_observation_data;
else
% KBR observations
intersat_observation_data = KBR_intersat_observation_data;
end
intersat_range_residuals = resrange;
intersat_rangerate_residuals = resrangerate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intersatellite ranging Output structure array
intersat_pod.intersat_observation_data    = intersat_observation_data;
intersat_pod.intersat_range_residuals     = intersat_range_residuals;
intersat_pod.intersat_rangerate_residuals = intersat_rangerate_residuals;

intersat_pod.KBR_rms_rangerate_residuals  = KBR_rms_res_rangerate;
intersat_pod.KBR_rms_range_residuals      = KBR_rms_res_range;

intersat_pod.LRI_rms_rangerate_residuals  = LRI_rms_res_rangerate;
intersat_pod.LRI_rms_range_residuals      = LRI_rms_res_range;

% Overall Output structure array
grace_pod_struct.grace1_pod   = grace1_pod;
grace_pod_struct.grace2_pod   = grace2_pod;
grace_pod_struct.intersat_pod = intersat_pod;
grace_pod_struct.intersat_ranging_functionals = intersat_ranging_functionals;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write Intersatellite Observation residuals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Mission Directory: GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis
% [OUT_fname_object_mjd, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_GRACE1);
% [status, message, messageid] = rmdir(OUT_fname_mission_mjd,'s');
% [status, message, messageid] = mkdir(OUT_fname_mission_mjd);
% 
% if LRI_obs_analysis_01 == 1
% data_matrix = LRI_rangerate_residuals;
% data_functional = 'LRI range-rate residuals';
% reference_frame = 'ICRF';
% [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% 
% data_matrix = LRI_range_residuals;
% data_functional = 'LRI range residuals';
% reference_frame = 'ICRF';
% [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% end
% 
% data_matrix = KBR_rangerate_residuals;
% data_functional = 'KBR range-rate residuals';
% reference_frame = 'ICRF';
% [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% 
% data_matrix = KBR_range_residuals;
% data_functional = 'KBR range residuals';
% reference_frame = 'ICRF';
% [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write Inter-Satellite Ranging functionals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % intersat_ranging_functionals
% % intersat_range_functionals = intersat_ranging_functionals.range ;
% % intersat_rangerate_functionals = intersat_ranging_functionals.rangerate ;
% % 
% % data_matrix = intersat_rangerate_functionals;
% % data_functional = 'Intersatellite range-rate functionals';
% % reference_frame = 'ICRF';
% % [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% % [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% % % GNP format
% % [georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% % [status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_mission_mjd);
% % 
% % data_matrix = intersat_range_functionals;
% % data_functional = 'Intersatellite range functionals';
% % reference_frame = 'ICRF';
% % [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% % [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% % % GNP format
% % [georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% % [status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_mission_mjd);
% 
% intersat_ranging_functionals = intersat_ranging_functionals.rangingdata ;
% data_matrix = intersat_ranging_functionals;
% data_functional = 'Intersatellite ranging functionals';
% reference_frame = 'ICRF';
% [georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write Statistics
%     rms_extorb(1,1:3) = rms_orbitalGA;
%     rms_extorb(2,1:3) = rms_orbitalGB;
%     rms_obs(1,1:3)    = rms_obs_GA;
%     rms_obs(2,1:3)    = rms_obs_GB;
%     rms_kbr = zeros(1,3);
%     rms_lri = zeros(1,3);
%     rms_kbr(1,2) = KBR_rms_resrangerate;
% if LRI_obs_analysis_01 == 1
%     rms_lri(1,2) = LRI_rms_resrangerate;
% end
% [georb_data_name] = write_georb_statistics(orbit_config_G1, orbit_config_G2, rms_obs, rms_extorb, rms_kbr, rms_lri);
% [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [OUT_fname_G1, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_GRACE1);
% [OUT_fname_G2, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G2,orbit_model_GRACE2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Move GRACE/GRACE-FO mission data files to one directory
%     [status,message,messageid] = movefile('*.out',OUT_fname_mission_mjd);
%     [status,message,messageid] = movefile('*.orb',OUT_fname_mission_mjd);
%     [status,message,messageid] = movefile(OUT_fname_G1,OUT_fname_mission_mjd);
%     [status,message,messageid] = movefile(OUT_fname_G2,OUT_fname_mission_mjd);
%     out_dir_name = OUT_fname_mission_mjd;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
