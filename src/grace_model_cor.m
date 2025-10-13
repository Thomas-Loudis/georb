function [grace_pod_struct, orbit_model_G1, orbit_model_G2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, gravity_sol_01, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: grace_model_cor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_model_struct  : Orbit model structure array 
% - ic_data_struct      : Initial Conditions structure array 
% - Xmatrix_flag (ic_apriori_01) :
%    0 : Aposteriori/Apriori values of IC Parameters are provided by the variable Xmatrix
%    1 : Corrections to the IC Parameters are provided by the variable Xmatrix
%   -1 : IC is not update; Xmatrix variables are ignored 
%
% Output arguments:
% - out_dir_name        : Output Data folder name 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                              26 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Structure arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orbit_model_G1  = grace_pod_struct.grace1_pod.orbit_model;
orbit_model_G2  = grace_pod_struct.grace2_pod.orbit_model;

orbit_config_G1 = grace_pod_struct.grace1_pod.orbit_config;
orbit_config_G2 = grace_pod_struct.grace2_pod.orbit_config;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter update
grav_param_01 = gravity_sol_01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer Calibration parameter update
acc_cal_update = 0;
% if grav_param_01 == 1 
% acc_cal_update = 1; 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions change
Xmatrix = Xmatrix_orbit1;

% GRACE-1/3 orbit model/configuration array update 
[orbit_config_G1,orbit_model_G1] = config_struct_update(orbit_config_G1, orbit_pod_mode, Xmatrix, Xmatrix_flag, acc_cal_update, grav_param_01,orbit_model_G1);

% Initial Conditions change
Xmatrix = Xmatrix_orbit2;

% GRACE-2/4 orbit model/configuration array update 
[orbit_config_G2,orbit_model_G2] = config_struct_update(orbit_config_G2, orbit_pod_mode, Xmatrix, Xmatrix_flag, acc_cal_update, grav_param_01,orbit_model_G2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orbit_model_G1.orbit_config = orbit_config_G1;
orbit_model_G2.orbit_config = orbit_config_G2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_pod strucutre array update 
grace_pod_struct.grace1_pod.orbit_model = orbit_model_G1;
grace_pod_struct.grace2_pod.orbit_model = orbit_model_G2;

grace_pod_struct.grace1_pod.orbit_config = orbit_config_G1;
grace_pod_struct.grace2_pod.orbit_config = orbit_config_G2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field solution update is handled individually via function gravity_param_aposteriori 
% [Xaposteriori, orbit_model_struct] = gravity_param_aposteriori(Xmatrix, ic_apriori_01, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

