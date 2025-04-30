function [OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_mission_grace
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
% Dr. Thomas Loudis Papanikolaou                            23 August  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 26/04/2025, Thomas Loudis Papanikolaou
%             Code extracted from function orbit_mission_grace and modfied
%             accordingly  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Structure arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Models (Force modelling matrices, Methods parameters, Data preprocessing)
orbit_model_matrix_GRACE1   = grace_pod_struct.grace1_pod.orbit_model;
orbit_model_matrix_GRACE2   = grace_pod_struct.grace2_pod.orbit_model;

% Orbit configuration (Input configuration variables)
orbit_config_G1             = grace_pod_struct.grace1_pod.orbit_config;
orbit_config_G2             = grace_pod_struct.grace2_pod.orbit_config;

% Orbit resutls (orbits, residuals, statistics)
orbit_matrices_G1           = grace_pod_struct.grace1_pod.orbit_matrices;
orbit_matrices_G2           = grace_pod_struct.grace2_pod.orbit_matrices;

% Inter-Satellite Ranging data matrices
intersat_pod = grace_pod_struct.intersat_pod;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit matrices in ICRF 
orbcGA                      = grace_pod_struct.grace1_pod.orbit_matrices.orbit_crf;
orbcGB                      = grace_pod_struct.grace2_pod.orbit_matrices.orbit_crf;

% Orbital frame residuals
rms_orbitalGA = grace_pod_struct.grace1_pod.orbit_matrices.rms_orbital;
rms_orbitalGB = grace_pod_struct.grace2_pod.orbit_matrices.rms_orbital;

% Observation matrix :: Pseudo-Observation GPS-based
sat1_OBS_matrix = grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix;
sat2_OBS_matrix = grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix;
% orbit_matrices.observation_residuals = OBS_residuals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE Inter-Satellite Data matrices
% intersat_observation_data = grace_pod_struct.intesat_pod_struct.intersat_observation_data;
intersat_observation_data       = grace_pod_struct.intersat_pod.intersat_observation_data ;
intersat_range_residuals        = grace_pod_struct.intersat_pod.intersat_range_residuals;
intersat_rangerate_residuals    = grace_pod_struct.intersat_pod.intersat_rangerate_residuals;

KBR_rms_resrangerate           = intersat_pod.KBR_rms_rangerate_residuals;
KBR_rms_res_range              = intersat_pod.KBR_rms_range_residuals ;

LRI_rms_resrangerate           = intersat_pod.LRI_rms_rangerate_residuals;
LRI_rms_res_range              = intersat_pod.LRI_rms_range_residuals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRI_rangerate_residuals         = intersat_rangerate_residuals;
LRI_range_residuals             = intersat_range_residuals;
KBR_rangerate_residuals         = LRI_rangerate_residuals;
KBR_range_residuals             = LRI_range_residuals;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite Observations 
LRI_obs_analysis_01 = 0;
param_keyword = 'LRI';
[test_intersat_obs_flag] = strcmp(intersat_obs_flag,param_keyword);
if test_intersat_obs_flag == 1
    LRI_obs_analysis_01 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orbit_config_struct_GRACE1 = orbit_config_G1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Orbit Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE1
orbc = orbcGA;
orbce = sat1_OBS_matrix;
orbc_pos = orbc(:,1:4);
% Celestial Reference Frame (GCRS)
[dorbc,rms_orbc,orbc_common] = compstat(orbce,orbc_pos);
rms_obs_GA = rms_orbc;
OBS_residuals = dorbc;

% Update central structure arrays
orbit_matrices_G1.observation_residuals = OBS_residuals;
grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals = OBS_residuals;

% GRACE2
orbc = orbcGB;
orbce = sat2_OBS_matrix;
orbc_pos = orbc(:,1:4);
% Celestial Reference Frame (GCRS)
[dorbc,rms_orbc,orbc_common] = compstat(orbce,orbc_pos);
rms_obs_GB = rms_orbc;
OBS_residuals = dorbc;

% Update central structure arrays
orbit_matrices_G2.observation_residuals = OBS_residuals;
grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals = OBS_residuals;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write orbit data in georb format and move to folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write orbit data of GRACE1
[OUT_data_foldername_G1] = write_georb_data_products(orbit_config_G1, orbit_model_matrix_GRACE1, orbit_matrices_G1); 

% Write orbit data of GRACE1
[OUT_data_foldername_G2] = write_georb_data_products(orbit_config_G2, orbit_model_matrix_GRACE2, orbit_matrices_G2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Intersatellite Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mission Directory: GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis
[OUT_fname_G1, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_matrix_GRACE1);
[OUT_fname_G2, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G2,orbit_model_matrix_GRACE2);
[status, message, messageid] = rmdir(OUT_fname_mission_mjd,'s');
[status, message, messageid] = mkdir(OUT_fname_mission_mjd);

if LRI_obs_analysis_01 == 1
data_matrix = LRI_rangerate_residuals;
data_functional = 'LRI range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = LRI_range_residuals;
data_functional = 'LRI range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
end

data_matrix = KBR_rangerate_residuals;
data_functional = 'KBR range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = KBR_range_residuals;
data_functional = 'KBR range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Statistics
    rms_extorb(1,1:3) = rms_orbitalGA;
    rms_extorb(2,1:3) = rms_orbitalGB;
    rms_obs(1,1:3)    = rms_obs_GA;
    rms_obs(2,1:3)    = rms_obs_GB;
    rms_kbr = zeros(1,3);
    rms_lri = zeros(1,3);
    rms_kbr(1,2) = KBR_rms_resrangerate;
if LRI_obs_analysis_01 == 1
    rms_lri(1,2) = LRI_rms_resrangerate;
end
    [georb_data_name] = write_georb_statistics(orbit_config_G1, orbit_config_G2, rms_obs, rms_extorb, rms_kbr, rms_lri);
    [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations matrices
NEQ_N = grace_pod_struct.normal_equations_Nmatrix;
NEQ_u = grace_pod_struct.normal_equations_umatrix;
NEQ_N_reduced = grace_pod_struct.normal_equations_Nmatrix_reduced;
NEQ_u_reduced = grace_pod_struct.normal_equations_umatrix_reduced;

data_matrix = NEQ_N;
data_functional = 'NEQ_N_matrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = NEQ_u;
data_functional = 'NEQ_u_matrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

if 1 < 0
% NEQ_N_reduced, NEQ_u_reduced
data_matrix = NEQ_N_reduced;
data_functional = 'NEQ_N_reduced_matrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = NEQ_u_reduced;
data_functional = 'NEQ_u_reduced_matrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move GRACE/GRACE-FO mission data files to one directory
[status,message,messageid] = movefile('*.out',OUT_fname_mission_mjd);
[status,message,messageid] = movefile('*.orb',OUT_fname_mission_mjd);
% [status,message,messageid] = movefile(OUT_fname_G1,OUT_fname_mission_mjd);
% [status,message,messageid] = movefile(OUT_fname_G2,OUT_fname_mission_mjd);
[status,message,messageid] = movefile(OUT_data_foldername_G1,OUT_fname_mission_mjd);
[status,message,messageid] = movefile(OUT_data_foldername_G2,OUT_fname_mission_mjd);
out_dir_name = OUT_fname_mission_mjd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
