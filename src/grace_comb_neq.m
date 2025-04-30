function [NEQ_N, NEQ_u, NEQ_N_reduced, NEQ_u_reduced, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xcommon_NEQreduced, grace_pod_struct] = grace_comb_neq(grace_pod_struct)

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
% 22/04/2025, Thomas Loudis Papanikolaou
%             Code extracted from function orbit_mission_grace and modified accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data from grace_pod structure array
orbit_model_GRACE1 = grace_pod_struct.grace1_pod.orbit_model;
orbit_model_GRACE2 = grace_pod_struct.grace2_pod.orbit_model;
orbit_model_struct = orbit_model_GRACE1;

% orbit_config_struct_GRACE1 = grace_pod_struct.grace1_pod.orbit_config;
% orbit_config_struct_GRACE2 = grace_pod_struct.grace2_pod.orbit_config;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Parameter estimation settings (TEMP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observation_model_matrix = orbit_model_struct.observation_model_matrix;
gravity_field_determination_matrix = orbit_model_struct.gravity_field_determination_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMBESTIM_combparamestim_01 = observation_model_matrix.COMBESTIM_combparamestim_01;

% Number of estimation occurences (If 0 :: Combined-Estimator not applied)   
COMBESTIM_Nestim_comb = observation_model_matrix.COMBESTIM_Nestim_comb;

% Combined solution Observations:  
COMBESTIM_OBS = observation_model_matrix.COMBESTIM_OBS;

% Intersatellite Observations :: Laser (LRI) or K-band (KBR)
COMBESTIM_intersat_obs = observation_model_matrix.COMBESTIM_intersat_obs;

% Weighted Estimation solution approaches
COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data from grace_pod structure array
orbcGA                      = grace_pod_struct.grace1_pod.orbit_matrices.orbit_crf;
orbcGB                      = grace_pod_struct.grace2_pod.orbit_matrices.orbit_crf;
sat1_OBS_matrix             = grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix;
sat2_OBS_matrix             = grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix;

sat1_OBS_residuals          = grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals;
sat2_OBS_residuals          = grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals;

sat1_veqZ_matrix             = grace_pod_struct.grace1_pod.orbit_matrices.state_transition_matrix;
sat1_veqP_matrix             = grace_pod_struct.grace1_pod.orbit_matrices.sensitivity_matrix;

sat2_veqZ_matrix             = grace_pod_struct.grace2_pod.orbit_matrices.state_transition_matrix;
sat2_veqP_matrix             = grace_pod_struct.grace2_pod.orbit_matrices.sensitivity_matrix;

intersat_observation_data    = grace_pod_struct.intersat_pod.intersat_observation_data;
intersat_range_residuals     = grace_pod_struct.intersat_pod.intersat_range_residuals;
intersat_rangerate_residuals = grace_pod_struct.intersat_pod.intersat_rangerate_residuals;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations Weights 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted estimation solution approaches
% observation_model_matrix = orbit_model_GRACE1.observation_model_matrix;
% COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight;
weight_sol_opt = COMBESTIM_weight;
    
% Weights based on errors / residuals from step-1 orbit parameter estimation  
Cv_sat1_obs = [sat1_OBS_residuals(:,2); sat1_OBS_residuals(:,3); sat1_OBS_residuals(:,4)];
Cv_sat2_obs = [sat2_OBS_residuals(:,2); sat2_OBS_residuals(:,3); sat2_OBS_residuals(:,4)];  
% LRI range and range-rate obseravations
Cv_LRI_range = intersat_range_residuals(:,2);
Cv_LRI_rangerate = intersat_rangerate_residuals(:,2);
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];

% Fixed weights to sigma defined values
if weight_sol_opt == 3
sigma_obsorb = 5 * 10^-2;
sigma_range = 5 * 10^-3;
sigma_rangerate = 1 * 10^-6;
    
[d1_obs, d2_obs] = size(Cv_sat1_obs);
[d1_obs2, d2_obs2] = size(Cv_sat2_obs);
[d1_lri, d2_lri] = size(Cv_LRI_range);

Cv_sat1_obs = sigma_obsorb + zeros(d1_obs,1);
Cv_sat2_obs = sigma_obsorb + zeros(d1_obs2,1);
Cv_LRI_range = sigma_range + zeros(d1_lri,1);
Cv_LRI_rangerate = sigma_rangerate + zeros(d1_lri,1);
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];

elseif weight_sol_opt == 1
% Identity matrix
Cv_sat1_obs = 1;
Cv_sat2_obs = 1;
Cv_LRI_range = 1;
Cv_LRI_rangerate = 1;
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-4 matrices for combined estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite observatios w.r.t. GRACE-4
[Xmatrix_LRI_sat2, Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = estimator_orbit_intersat(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, Cv_LRI);

% Orbit pseudo-observations
[Xmatrix_obsorb_grace2,Xmatrix_alt_obsorb_grace2,Wmatrix_obsorb_grace2,Amatrix_obsorb_grace2, Cx, Cv, NEQn_grace2, NEQu_grace2] = estimator_orbit (orbcGB, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, Cv_sat2_obs,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-3 matrices for combined estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite observatios w.r.t. GRACE-4
[Xmatrix_LRI_sat1, Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = estimator_orbit_intersat(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, Cv_LRI);

% Orbit pseudo-observations
[Xmatrix_obsorb_grace1,Xmatrix_alt_obsorb_grace1,Wmatrix_obsorb_grace1,Amatrix_obsorb_grace1, Cx, Cv, NEQn_grace1, NEQu_grace1] = estimator_orbit (orbcGA, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, Cv_sat1_obs,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gravity_field_determination_matrix = orbit_model_GRACE1.gravity_field_determination_matrix;
COMBESTIM_gravity_step = gravity_field_determination_matrix.COMBESTIM_gravity_step;
gfm_struct = orbit_model_GRACE1.gravity_field;
Nparam_common = 0;
if COMBESTIM_gravity_step > 0
grav_paramestim_yn = gfm_struct.param_estim_yn;
test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');   
if test_grav_paramestim_yn == 1   
N_param_GRAV = gfm_struct.parameters_number;
Nparam_common = N_param_GRAV;
end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions of design matrix and Number of parameters, observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d1 d2]= size(Amatrix_obsorb_grace1);
[d3 d4]= size(Amatrix_obsorb_grace2);
Nparam_grace1 = d2;
Nparam_grace2 = d4;
Nobs_sat1 = d1;
Nobs_sat2 = d3;
Nparam_COMB = Nparam_grace1 + Nparam_grace2 - Nparam_common;

[d5 d6]= size(Amatrix_rangerate_sat1);
[d7 d8]= size(Amatrix_rangerate_sat2);
Nobs_rangerate = d7;
Nparam_rangerate = d8;

[ d9 d10]= size(Amatrix_range_sat1);
[d11 d12]= size(Amatrix_range_sat2);
Nobs_range = d9;
Nparam_range = d10;

N_obs_comb = Nobs_sat1 + Nobs_sat2 + Nobs_rangerate + Nobs_range;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amatrix = zeros(N_obs_comb , Nparam_COMB);
col_orbit1_1 = 1;
col_orbit1_2 = Nparam_grace1-Nparam_common;
col_orbit2_1 = col_orbit1_2 + 1;
col_orbit2_2 = col_orbit1_2 + Nparam_grace2-Nparam_common;
col_common_1 = col_orbit2_2 + 1;
col_common_2 = col_orbit2_2 + Nparam_common;

% Selected parameters to be estimated :: All or Golbal parameters (common)
param_col_start = 1;
param_col_end = Nparam_COMB;

% Orbit 1
A_orbit1 = zeros(Nobs_sat1 , param_col_end - (param_col_start-1));
% A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1; 
A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1(:,1: col_orbit1_2);
if Nparam_common > 0
A_orbit1(1 : Nobs_sat1 , col_common_1 : col_common_2) = Amatrix_obsorb_grace1(:, Nparam_grace1-Nparam_common+1 : Nparam_grace1-Nparam_common+Nparam_common);
end

% Orbit 2
A_orbit2 = zeros(Nobs_sat2 , param_col_end - (param_col_start-1));
% A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2; 
A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2(:,1:Nparam_grace2-Nparam_common);
if Nparam_common > 0
A_orbit2(1 : Nobs_sat2 , col_common_1 : col_common_2) = Amatrix_obsorb_grace2(:,Nparam_grace2-Nparam_common+1 : Nparam_grace2-Nparam_common+Nparam_common);
end

% Range-rate
col_common_0 = Nparam_rangerate-Nparam_common;
A_rangerate(1 : Nobs_rangerate , col_orbit1_1 : col_orbit1_2) = Amatrix_rangerate_sat1(:,1:col_common_0);
A_rangerate(1 : Nobs_rangerate , col_orbit2_1 : col_orbit2_2) = Amatrix_rangerate_sat2(:,1:col_common_0);
if Nparam_common > 0
A_rangerate(1 : Nobs_rangerate , col_common_1 : col_common_2) = Amatrix_rangerate_sat1(:,col_common_0+1:col_common_0+Nparam_common) + Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_common); 
end

% Range
col_common_0 = Nparam_range-Nparam_common;
A_range(1 : Nobs_range, col_orbit1_1 : col_orbit1_2) = Amatrix_range_sat1(:,1:col_common_0);
A_range(1 : Nobs_range, col_orbit2_1 : col_orbit2_2) = Amatrix_range_sat2(:,1:col_common_0);
if Nparam_common > 0
A_range(1 : Nobs_range, col_common_1 : col_common_2) = Amatrix_range_sat1(:,col_common_0+1:col_common_0+Nparam_common) + Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_common);
end

% Reduced-Observations matrix
b_orbit1 = Wmatrix_obsorb_grace1;
b_orbit2 = Wmatrix_obsorb_grace2;
b_rangerate = Wmatrix_rangerate_sat2;
b_range = Wmatrix_range_sat2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Approach-3
% Individual NEQ matrices per block/type of observations
% Orbit 1
[Xestim_orbit1, NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq_sol(A_orbit1, b_orbit1, Cv_sat1_obs);
% Orbit 2
[Xestim_orbit2, NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq_sol(A_orbit2, b_orbit2, Cv_sat2_obs);
% Range-rate
[Xestim_rangerate, NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq_sol(A_rangerate, b_rangerate, Cv_LRI_rangerate);
% Range
[Xestim_range, NEQ_N_range, NEQ_u_range] = estimator_neq_sol(A_range, b_range, Cv_LRI_range);

% Combined solution Observations:  
% Normal Equations :: Sum matrices
if COMBESTIM_OBS == 1
NEQ_N = NEQ_N_orbit1;
NEQ_u = NEQ_u_orbit1;

elseif COMBESTIM_OBS == 2
NEQ_N = NEQ_N_orbit2;
NEQ_u = NEQ_u_orbit2;
    
elseif COMBESTIM_OBS == 3
NEQ_N = NEQ_N_orbit1 + NEQ_N_orbit2;
NEQ_u = NEQ_u_orbit1 + NEQ_u_orbit2;
    
elseif COMBESTIM_OBS == 4
NEQ_N = NEQ_N_range;
NEQ_u = NEQ_u_range;
  
elseif COMBESTIM_OBS == 5
NEQ_N = NEQ_N_rangerate;
NEQ_u = NEQ_u_rangerate;

elseif COMBESTIM_OBS == 6
NEQ_N = NEQ_N_rangerate + NEQ_N_range;
NEQ_u = NEQ_u_rangerate + NEQ_u_range;

elseif COMBESTIM_OBS == 7 
NEQ_N = NEQ_N_orbit1 + NEQ_N_orbit2 + NEQ_N_range;
NEQ_u = NEQ_u_orbit1 + NEQ_u_orbit2 + NEQ_u_range;

elseif COMBESTIM_OBS == 8
NEQ_N = NEQ_N_orbit1 + NEQ_N_orbit2 + NEQ_N_rangerate + NEQ_N_range;
NEQ_u = NEQ_u_orbit1 + NEQ_u_orbit2 + NEQ_u_rangerate + NEQ_u_range;

elseif COMBESTIM_OBS == 9
NEQ_N = NEQ_N_orbit1 + NEQ_N_orbit2 + NEQ_N_rangerate;
NEQ_u = NEQ_u_orbit1 + NEQ_u_orbit2 + NEQ_u_rangerate;
end

% Weighted Least Squares method solution
tol2 = 30;
Xmatrix3 = lsqminnorm(NEQ_N, NEQ_u, tol2);
Xmatrix = Xmatrix3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced NEQ matrices :: pre-elimination of orbit arc-related parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_orbparam  = (Nparam_grace1 - Nparam_common) + (Nparam_grace2 - Nparam_common);
n_gravparam = Nparam_common;
[Xmatrix_NEQ_reduced, NEQ_N_reduced, NEQ_u_reduced] = neq_reduced(NEQ_N, NEQ_u, n_orbparam, n_gravparam);
Xmatrix_GRAV_NEQ_reduced = Xmatrix_NEQ_reduced;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit and Gravity Field parameters estimated corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Parameters estimated corrections    
Xmatrix_orbit1 = Xmatrix(1 : col_orbit1_2 , 1);
Xmatrix_orbit2 = Xmatrix(col_orbit2_1 : col_orbit2_2 , 1);   
% Gravity parameters estimated corrections    
Xmatrix_gravparam = Xmatrix(col_common_1 : col_common_2 , 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xcommon = Xmatrix_gravparam;
Xcommon_NEQreduced = Xmatrix_NEQ_reduced;

% Structure array update :: param_aposteriori
% sat1_Xparam_aposteriori     = grace_pod_struct.grace1_pod.orbit_matrices.param_aposteriori;
grace_pod_struct.grace1_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit1;
grace_pod_struct.grace2_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit2;
