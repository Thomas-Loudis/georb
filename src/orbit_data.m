function [accelerometer_struct, observation_struct, KBR_intersat_struct, LRI_intersat_struct] = orbit_data (cfg_fname, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit data : Read the data required for the orbit modelling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  26 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE mission: Accelerometry Data (Calibration parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'acc_data';
[acc_data] = read_param_cfg(cfg_fname,param_keyword);
test_accelerometer = strcmp(acc_data,'y');
if test_accelerometer == 1        
% GRACE Accelerometer data preprocessing
[acc_cal_param, acc_dpint, acc1b_array, sca1b_array, sca_dpint, acc_cal_bias_yn, acc_scale_type, accelerometer_struct] = prm_grace_data(cfg_fname); 
% Number of accelerometer calibration parameters    
[d1, d2] = size(acc_cal_param);
Nparam_ACC_CAL = d1;
else
    % Accelerometer data use y/n 
    accelerometer_struct.effect_yn = acc_data;
    % Accelerometer calibration parameters estimation y/n 
    accelerometer_struct.param_estim_yn = 'n'; 
    Nparam_ACC_CAL = 0;
    accelerometer_struct.ACC1B_data_array = 0;
    accelerometer_struct.SCA1B_data_array = 0;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation data :: Pseudo-Observations to orbits 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'orbit_mode';
[orbit_mode] = read_param_cfg(cfg_fname,param_keyword);
test_pod_mode = strcmp(orbit_mode,'orbit_determination');
if test_pod_mode == 1   
% Pseudo-Observations via GPS-based kinematic orbit data
[obsorbc,obsorbt,obsorbc_ext,obsorbt_ext,obsorbc_full,obsorbt_full,COVobs,COVPform, observation_struct] = orbit_obs(cfg_fname, orbit_model_struct);
else    
observation_struct.obs_orbit_crf = zeros(4,4);
observation_struct.obs_orbit_trf = zeros(4,4);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation data :: Intersatellite Ranging data :: KBR & LRI data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB mode 
param_keyword = 'georb_mode';
% [georb_mode] = read_param_file(main_config_fname,param_keyword);
[georb_mode] = read_param_cfg(cfg_fname,param_keyword);

% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_objects_mission';
% [orbiting_object_name] = read_param_file(main_config_fname,param_keyword);
[orbiting_object_name] = read_param_cfg(cfg_fname,param_keyword);

% Mode: georb_mode :: 'orbit_mission' 
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Satellite missions cases :: 
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission');

if test_orbit_mission == 1 && ( test_mission_grace == 1 || test_mission_gracefo == 1 )
% KBR data reading 
intersat_type = 'KBR_data';
[KBR_intersat_struct] = grace_data_intersatellite (cfg_fname, intersat_type);

% LRI data read
if test_mission_gracefo == 1
intersat_type = 'LRI_data';
[LRI_intersat_struct] = grace_data_intersatellite (cfg_fname, intersat_type);
else
   LRI_intersat_struct.range = 0;
   LRI_intersat_struct.rangerate = 0;
   LRI_intersat_struct.rangeacceleration = 0;    
end

else
   KBR_intersat_struct.range = 0;
   KBR_intersat_struct.rangerate = 0;
   KBR_intersat_struct.rangeacceleration = 0;

   LRI_intersat_struct.range = 0;
   LRI_intersat_struct.rangerate = 0;
   LRI_intersat_struct.rangeacceleration = 0;       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

