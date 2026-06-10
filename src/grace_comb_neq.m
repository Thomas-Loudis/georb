function [NEQ_N, NEQ_u, NEQ_N_reduced, NEQ_u_reduced, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xcommon, Xcommon_NEQreduced, grace_pod_struct] = grace_comb_neq(grace_pod_struct, design_matrix_scale_orbit, design_matrix_param_orbit, NEQ_write, rangerate_empirical_struct) 

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
 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Combined Parameter estimation settings (TEMP) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% observation_model_matrix = orbit_model_struct.observation_model_matrix; 
% gravity_field_determination_matrix = orbit_model_struct.gravity_field_determination_matrix; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% COMBESTIM_combparamestim_01 = observation_model_matrix.COMBESTIM_combparamestim_01; 
%  
% % Number of estimation occurences (If 0 :: Combined-Estimator not applied)    
% COMBESTIM_Nestim_comb = observation_model_matrix.COMBESTIM_Nestim_comb; 
%  
% % Combined solution Observations:   
% COMBESTIM_OBS = observation_model_matrix.COMBESTIM_OBS; 
%  
% % Intersatellite Observations :: Laser (LRI) or K-band (KBR) 
% COMBESTIM_intersat_obs = observation_model_matrix.COMBESTIM_intersat_obs; 
%  
% % Weighted Estimation solution approaches 
% COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Combined solution Observations:   
COMBESTIM_OBS = orbit_model_struct.observation_model_matrix.COMBESTIM_OBS; 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Read data from grace_pod structure array 
orbcGA                       = grace_pod_struct.grace1_pod.orbit_matrices.orbit_crf; 
orbcGB                       = grace_pod_struct.grace2_pod.orbit_matrices.orbit_crf; 
sat1_OBS_matrix              = grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix; 
sat2_OBS_matrix              = grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix; 
 
sat1_OBS_residuals           = grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals; 
sat2_OBS_residuals           = grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals; 
 
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
[sigma_obs_orbit1, sigma_obs_orbit2, sigma_intersat_range, sigma_intersat_rangerate] = intersat_weight_function(grace_pod_struct); 
Cv_sat1_obs = sigma_obs_orbit1; 
Cv_sat2_obs = sigma_obs_orbit2; 
Cv_LRI_range = sigma_intersat_range; 
Cv_LRI_rangerate = sigma_intersat_rangerate; 
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GRACE-4 matrices for combined estimator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inter-satellite observations w.r.t. GRACE-2/-4 
% [Xmatrix_LRI_sat2, Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = estimator_orbit_intersat(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, Cv_LRI); 
% [Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = neq_intersat_ranging(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, Cv_LRI);  
[Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = neq_intersat_ranging(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, Cv_LRI, rangerate_empirical_struct);  
 
% Orbit pseudo-observations 
% [Xmatrix_obsorb_grace2,Xmatrix_alt_obsorb_grace2,Wmatrix_obsorb_grace2,Amatrix_obsorb_grace2, Cx, Cv, NEQn_grace2, NEQu_grace2] = estimator_orbit (orbcGB, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, Cv_sat2_obs,1); 
[NEQn_grace2, NEQu_grace2, Amatrix_obsorb_grace2, Wmatrix_obsorb_grace2] = neq_orbit (orbcGB, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, Cv_sat2_obs); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GRACE-3 matrices for combined estimator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inter-satellite observations w.r.t. GRACE-1/-3 
% [Xmatrix_LRI_sat1, Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = estimator_orbit_intersat(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, Cv_LRI); 
% [Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = neq_intersat_ranging(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, Cv_LRI);  
[Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = neq_intersat_ranging(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, Cv_LRI, rangerate_empirical_struct);  
 
% Orbit pseudo-observations 
% [Xmatrix_obsorb_grace1,Xmatrix_alt_obsorb_grace1,Wmatrix_obsorb_grace1,Amatrix_obsorb_grace1, Cx, Cv, NEQn_grace1, NEQu_grace1] = estimator_orbit (orbcGA, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, Cv_sat1_obs,1); 
[NEQn_grace1, NEQu_grace1, Amatrix_obsorb_grace1, Wmatrix_obsorb_grace1] = neq_orbit (orbcGA, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, Cv_sat1_obs); 
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
% Range-Rate Empirical parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% rangerate_emp_01 = 1 
rangerate_emp_01       = rangerate_empirical_struct.rangerate_empirical_function_01; 
rangerate_emp_param    = rangerate_empirical_struct.rangerate_empirical_parameters; 
Nparam_rangerate_cycle = rangerate_empirical_struct.Nparam_rangerate_cycle;   
Nparam_rangerate_all   = rangerate_empirical_struct.Nparam_rangerate      ; 
 
% if rangerate_emp_01 == 1 
% [d1, d2] = size(Amatrix_obsorb_grace1); 
% [d3, d4] = size(Amatrix_obsorb_grace2); 
% [d5, d6] = size(rangerate_emp_param); 
%  
% Aq1 = zeros(d1,d5); 
% Aq2 = zeros(d3,d5); 
%  
% %    [rangerate_emp_obs, PD_rangerate_emp] = intersat_rangerate_emp(mjd, mjd0, zA, zB, rangerate_emp_param); 
% %    PD_range_emp = PD_rangerate_emp - PD_rangerate_emp; 
%  
% Amatrix_obsorb_grace1 = [Amatrix_obsorb_grace1 Aq1]; 
% Amatrix_obsorb_grace2 = [Amatrix_obsorb_grace2 Aq2]; 
% end 
 
if rangerate_emp_01 == 1 
% [Nparam_obs_rangerate, d6] = size(rangerate_emp_param); 
Nparam_obs_rangerate = Nparam_rangerate_all; 
else 
Nparam_obs_rangerate = 0;     
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
% Overall parameters number of the combined parameter estimation including gravity 
% parameters (global parameters) 
Nparam_COMB = Nparam_grace1 + Nparam_grace2 - Nparam_common + Nparam_obs_rangerate; 
% Number of orbital parameters (individual arc parameters) 
Nparam_orbital = Nparam_COMB - Nparam_common; 
 
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
% Design matrices' Collumns index numbers 
col_orbit1_1 = 1; 
col_orbit1_2 = Nparam_grace1-Nparam_common; 
 
col_orbit2_1 = col_orbit1_2 + 1; 
col_orbit2_2 = col_orbit1_2 + Nparam_grace2-Nparam_common; 
 
col_obs_rangerate_1 = col_orbit2_2 + Nparam_obs_rangerate / Nparam_obs_rangerate; 
col_obs_rangerate_2 = col_orbit2_2 + Nparam_obs_rangerate; 
 
col_common_1 = col_orbit2_2 + Nparam_obs_rangerate + 1; 
col_common_2 = col_orbit2_2 + Nparam_obs_rangerate + Nparam_common; 
 
% Selected parameters to be estimated :: All or Global parameters (common) 
param_col_start = 1; 
param_col_end = Nparam_COMB; 
 
% Orbit 1 
A_orbit1 = zeros(Nobs_sat1 , param_col_end - (param_col_start-1)); 
% Orbital Elements 
A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1(:,1: col_orbit1_2); 
% Common parameters :: Gravity 
if Nparam_common > 0 
A_orbit1(1 : Nobs_sat1 , col_common_1 : col_common_2) = Amatrix_obsorb_grace1(:, Nparam_grace1-Nparam_common+1 : Nparam_grace1-Nparam_common+Nparam_common); 
end 
 
% Orbit 2 
A_orbit2 = zeros(Nobs_sat2 , param_col_end - (param_col_start-1)); 
A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2(:,1:Nparam_grace2-Nparam_common); 
if Nparam_common > 0 
A_orbit2(1 : Nobs_sat2 , col_common_1 : col_common_2) = Amatrix_obsorb_grace2(:,Nparam_grace2-Nparam_common+1 : Nparam_grace2-Nparam_common+Nparam_common); 
end 
 
% Range-rate 
col_common_0 = Nparam_rangerate - Nparam_common - Nparam_obs_rangerate; 
A_rangerate(1 : Nobs_rangerate , col_orbit1_1 : col_orbit1_2) = Amatrix_rangerate_sat1(:,1:col_common_0); 
A_rangerate(1 : Nobs_rangerate , col_orbit2_1 : col_orbit2_2) = Amatrix_rangerate_sat2(:,1:col_common_0); 
% Observation Empirical Parameters 
if Nparam_obs_rangerate > 0 
A_rangerate(1 : Nobs_rangerate , col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_rangerate_sat1(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate); % + Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_obs_rangerate);  
end 
% Common parameters :: Gravity 
if Nparam_common > 0 
A_rangerate(1 : Nobs_rangerate , col_common_1 : col_common_2) = Amatrix_rangerate_sat1(:,col_common_0+1:col_common_0+Nparam_common) + Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_common);  
end 
 
% Range 
col_common_0 = Nparam_range - Nparam_common - Nparam_obs_rangerate; 
A_range(1 : Nobs_range, col_orbit1_1 : col_orbit1_2) = Amatrix_range_sat1(:,1:col_common_0); 
A_range(1 : Nobs_range, col_orbit2_1 : col_orbit2_2) = Amatrix_range_sat2(:,1:col_common_0); 
if Nparam_obs_rangerate > 0 
A_range(1 : Nobs_range, col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_range_sat1(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate); % + Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_obs_rangerate); 
end 
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
if NEQ_write == 1 
save A_orbit1.neq  A_orbit1 -ASCII -double 
save A_orbit2.neq  A_orbit2 -ASCII -double 
save A_rangerate.neq  A_rangerate -ASCII -double 
save A_range.neq  A_range -ASCII -double 
 
save b_orbit1.neq  b_orbit1 -ASCII -double 
save b_orbit2.neq  b_orbit2 -ASCII -double 
save b_rangerate.neq  b_rangerate -ASCII -double 
save b_range.neq  b_range -ASCII -double 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Scaling of Design Matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% design_matrix_scale_orbit 
if design_matrix_scale_orbit == 1 
scale_orbit_pos = 3 * 10^6; 
scale_orbit_vel = 3 * 10^3; 
% Orbit 1 
A_orbit1(:,col_orbit1_1 : col_orbit1_1 -1+3)         = A_orbit1(:,col_orbit1_1 : col_orbit1_1 -1+3) * scale_orbit_pos; 
A_orbit1(:,col_orbit1_1 -1+4 : col_orbit1_1 -1+6)    = A_orbit1(:,col_orbit1_1 -1+4 : col_orbit1_1 -1+6) * scale_orbit_vel; 
% Orbit 2 
A_orbit2(:,col_orbit2_1 : col_orbit2_1 -1+3)        = A_orbit2(:,col_orbit2_1 : col_orbit2_1 -1 + 3) * scale_orbit_pos; 
A_orbit2(:,col_orbit2_1 -1+4 : col_orbit2_1 -1+6)   = A_orbit2(:,col_orbit2_1 -1+4 : col_orbit2_1 -1+6) * scale_orbit_vel; 
 
% Range-Rate 
A_rangerate(:,col_orbit1_1 : col_orbit1_1 -1+3)        = A_rangerate(:,col_orbit1_1 : col_orbit1_1 -1+3)      * scale_orbit_pos; 
A_rangerate(:,col_orbit1_1-1 +4 : col_orbit1_1-1 +6)   = A_rangerate(:,col_orbit1_1 -1+4 : col_orbit1_1 -1+6) * scale_orbit_vel; 
A_rangerate(:,col_orbit2_1 : col_orbit2_1-1 +3)        = A_rangerate(:,col_orbit2_1 : col_orbit2_1 -1 + 3)    * scale_orbit_pos; 
A_rangerate(:,col_orbit2_1-1 +4 : col_orbit2_1-1 +6)   = A_rangerate(:,col_orbit2_1 -1+4 : col_orbit2_1 -1+6) * scale_orbit_vel; 
 
% Range 
A_range(:,col_orbit1_1 : col_orbit1_1 -1+3)        = A_range(:,col_orbit1_1 : col_orbit1_1 -1+3)      * scale_orbit_pos; 
A_range(:,col_orbit1_1-1 +4 : col_orbit1_1-1 +6)   = A_range(:,col_orbit1_1 -1+4 : col_orbit1_1 -1+6) * scale_orbit_vel; 
A_range(:,col_orbit2_1 : col_orbit2_1-1 +3)        = A_range(:,col_orbit2_1 : col_orbit2_1 -1 + 3)    * scale_orbit_pos; 
A_range(:,col_orbit2_1-1 +4 : col_orbit2_1-1 +6)   = A_range(:,col_orbit2_1 -1+4 : col_orbit2_1 -1+6) * scale_orbit_vel; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbital Parameters combination transformation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% design_matrix_param_orbit 
if design_matrix_param_orbit == 1 
% Orbital parameters of sate vector replacement     
% k1 = x_orbit1 - x_orbit2; 
% k2 = x_orbit1 + x_orbit2; 
 
% Orbit 1 
A_orbit1(:,col_orbit1_1 : col_orbit1_2) = A_orbit1(:,col_orbit1_1 : col_orbit1_2) * (1/2); 
A_orbit1(:,col_orbit2_1 : col_orbit2_2) = A_orbit1(:,col_orbit1_1 : col_orbit1_2); 
% Orbit 2 
A_orbit2(:,col_orbit2_1 : col_orbit2_2) = A_orbit2(:,col_orbit2_1 : col_orbit2_2) * (1/2); 
A_orbit2(:,col_orbit1_1 : col_orbit1_2) = A_orbit2(:,col_orbit2_1 : col_orbit2_2) * (-1); 
% Range-Rate 
% A_rangerate(:, col_orbit1_1 : col_orbit1_2) = A_rangerate(:, col_orbit1_1 : col_orbit1_2); 
A_rangerate(:, col_orbit2_1 : col_orbit2_2) = A_rangerate(:, col_orbit2_1 : col_orbit2_2) - A_rangerate(:, col_orbit2_1 : col_orbit2_2); % zeros; 
% Range 
% A_range(:, col_orbit1_1 : col_orbit1_2) = A_range(:, col_orbit1_1 : col_orbit1_2); 
A_range(:, col_orbit2_1 : col_orbit2_2) = A_range(:, col_orbit2_1 : col_orbit2_2) - A_range(:, col_orbit2_1 : col_orbit2_2); % zeros; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Normal Equations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. Approach-3 
% Individual NEQ matrices per block/type of observations 
% Orbit 1 
% [Xestim_orbit1, NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq_sol(A_orbit1, b_orbit1, Cv_sat1_obs); 
[NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq(A_orbit1, b_orbit1, Cv_sat1_obs); 
% Orbit 2 
% [Xestim_orbit2, NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq_sol(A_orbit2, b_orbit2, Cv_sat2_obs); 
[NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq(A_orbit2, b_orbit2, Cv_sat2_obs); 
% Range-rate 
% [Xestim_rangerate, NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq_sol(A_rangerate, b_rangerate, Cv_LRI_rangerate); 
[NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq(A_rangerate, b_rangerate, Cv_LRI_rangerate); 
% Range 
% [Xestim_range, NEQ_N_range, NEQ_u_range] = estimator_neq_sol(A_range, b_range, Cv_LRI_range); 
[NEQ_N_range, NEQ_u_range] = estimator_neq(A_range, b_range, Cv_LRI_range); 
 
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
 
[NEQn_d1 NEQn_d2] = size(NEQ_N); 
[NEQu_d1 NEQu_d2] = size(NEQ_u); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Weighted Least Squares solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
inv_id = 6; 
% if Nparam_obs_rangerate > 0 
% inv_id = 4 
% end 
[Xmatrix] = inv_ls(NEQ_N, NEQ_u, inv_id); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Covariance matrices 
sigma0 = 1; 
 
% Covariance matrix of parameters 
% inv_id = 6;  
[Nmatrix_inv] = inv_mat(NEQ_N, inv_id); 
Cx_matrix = sigma0^2 * Nmatrix_inv; 
Cx = diag(Cx_matrix); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Reduced NEQ matrices :: pre-elimination of orbit arc-related parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if Nparam_common > 0 
% n_orbparam  = (Nparam_grace1 - Nparam_common) + (Nparam_grace2 - Nparam_common) + Nparam_obs_rangerate; 
n_orbparam  = Nparam_orbital; 
n_gravparam = Nparam_common; 
[Xmatrix_NEQ_reduced, NEQ_N_reduced, NEQ_u_reduced] = neq_reduced(NEQ_N, NEQ_u, n_orbparam, n_gravparam); 
else 
Xmatrix_NEQ_reduced = 0; 
NEQ_N_reduced = 0; 
NEQ_u_reduced = 0; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Cond_no_prior_NEQn_grace1 = cond(NEQn_grace1) 
% Cond_no_prior_NEQn_grace2 = cond(NEQn_grace2) 
Cond_no_post_Nfinal = cond(NEQ_N); 
Cond_no_post_Nfinal_reduced = cond(NEQ_N_reduced); 
 
Condition_number_matrix(1,1) = Cond_no_post_Nfinal; 
Condition_number_matrix(2,1) = Cond_no_post_Nfinal_reduced; 
 
if NEQ_write == 1 
save NEQ_N_orbit1.neq  NEQ_N_orbit1 -ASCII -double 
save NEQ_u_orbit1.neq  NEQ_u_orbit1 -ASCII -double 
 
save NEQ_N_orbit2.neq  NEQ_N_orbit2 -ASCII -double 
save NEQ_u_orbit2.neq  NEQ_u_orbit2 -ASCII -double 
 
save NEQ_N_rangerate.neq  NEQ_N_rangerate -ASCII -double 
save NEQ_u_rangerate.neq  NEQ_u_rangerate -ASCII -double 
 
save NEQ_N_range.neq  NEQ_N_range -ASCII -double 
save NEQ_u_range.neq  NEQ_u_range -ASCII -double 
 
save NEQ_N_reduced.neq  NEQ_N_reduced -ASCII -double 
save NEQ_u_reduced.neq  NEQ_u_reduced -ASCII -double 
 
save NEQ_N.neq  NEQ_N -ASCII -double 
save NEQ_u.neq  NEQ_u -ASCII -double 
 
save COND_number_NEQ.neq  Condition_number_matrix -ASCII -double 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit and Gravity Field parameters estimated corrections 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit Parameters estimated corrections     
Xmatrix_orbit1 = Xmatrix(1 : col_orbit1_2 , 1); 
Xmatrix_orbit2 = Xmatrix(col_orbit2_1 : col_orbit2_2 , 1) ; 
 
if design_matrix_param_orbit == 1 
% Orbital parameters of state vector replacement back to original      
% k1 = x_orbit1 - x_orbit2; 
% k2 = x_orbit1 + x_orbit2; 
% x_orbit1 =  (k1 + k2) / 2;  
% x_orbit2 = (-k1 + k2) / 2; 
k1 = Xmatrix_orbit1; 
k2 = Xmatrix_orbit2; 
Xmatrix_orbit1 =  (k1 + k2) / 2 ; 
Xmatrix_orbit2 = (-k1 + k2) / 2 ; 
 
Xmatrix(1 : col_orbit1_2 , 1) = Xmatrix_orbit1; 
Xmatrix(col_orbit2_1 : col_orbit2_2 , 1) = Xmatrix_orbit2;   
end 
 
if design_matrix_scale_orbit == 1 
% Xmatrix_orbit1 = [scale_orbit_pos * Xmatrix(1 : col_orbit1_1 -1+3 , 1); ...  
%                   scale_orbit_vel * Xmatrix(col_orbit1_1 -1+4 : col_orbit1_1 -1+6 , 1) ]; 
% Xmatrix_orbit2 = [scale_orbit_pos * Xmatrix(col_orbit2_1 : col_orbit2_1 -1+3 , 1); ...  
%                   scale_orbit_vel * Xmatrix(col_orbit2_1 -1+4 : col_orbit2_1 -1+6 , 1) ] ; 
 
Xmatrix(1 : col_orbit1_1 -1+3 , 1)                  = scale_orbit_pos * Xmatrix(1 : col_orbit1_1 -1+3 , 1);  
Xmatrix(col_orbit1_1 -1+4 : col_orbit1_1 -1+6 , 1)  = scale_orbit_vel * Xmatrix(col_orbit1_1 -1+4 : col_orbit1_1 -1+6 , 1); 
 
Xmatrix(col_orbit2_1 : col_orbit2_1 -1+3 , 1)       = scale_orbit_pos * Xmatrix(col_orbit2_1 : col_orbit2_1 -1+3 , 1);  
Xmatrix(col_orbit2_1 -1+4 : col_orbit2_1 -1+6 , 1)  = scale_orbit_vel * Xmatrix(col_orbit2_1 -1+4 : col_orbit2_1 -1+6 , 1) ; 
 
Xmatrix_orbit1 = Xmatrix(1 : col_orbit1_2 , 1); 
Xmatrix_orbit2 = Xmatrix(col_orbit2_1 : col_orbit2_2 , 1) ; 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Gravity parameters estimated corrections     
if Nparam_common > 0 
Xmatrix_common = Xmatrix(col_common_1 : col_common_2 , 1); 
Xcommon             = Xmatrix_common; 
else 
Xcommon = 0; 
end 
Xcommon_NEQreduced  = Xmatrix_NEQ_reduced; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Structure array update :: param_aposteriori 
grace_pod_struct.grace1_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit1; 
grace_pod_struct.grace2_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit2; 
 
grace_pod_struct.grace1_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit1; 
grace_pod_struct.grace2_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit2; 
 
if NEQ_write == 1 
save Xcommon_grav.neq  Xcommon -ASCII -double 
save Xcommon_grav_NEQreduced.neq  Xcommon_NEQreduced -ASCII -double 
end 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range-Rate Empirical Parameters update 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if Nparam_obs_rangerate > 0 
% col_obs_rangerate_1 = col_orbit2_2 + Nparam_obs_rangerate / Nparam_obs_rangerate; 
% col_obs_rangerate_2 = col_orbit2_2 + Nparam_obs_rangerate; 
Xmatrix_obs_rangerate = Xmatrix(col_obs_rangerate_1 : col_obs_rangerate_2 , 1) 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range-Rate Empirical Parameters matrix update 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% rangerate_emp_param_apriori = rangerate_empirical_struct.rangerate_empirical_parameters; 
% rangerate_empirical_struct.rangerate_empirical_parameters  = Xmatrix_obs_rangerate + rangerate_emp_param_apriori; 
 
% rangerate_emp_01       = rangerate_empirical_struct.rangerate_empirical_function_01; 
Nparam_rangerate_cycle = rangerate_empirical_struct.Nparam_rangerate_cycle;   
Nparam_rangerate_all   = rangerate_empirical_struct.Nparam_rangerate      ; 
 
% Apriori values  
rangerate_emp_param_apriori_matrix = rangerate_empirical_struct.rangerate_empirical_parameters; 
Num_time_collumns = 2; 
 
% Aposteriori values 
rangerate_emp_param_aposteriori_matrix = rangerate_emp_param_apriori_matrix; 
 
[d1,d2] = size(rangerate_emp_param_apriori_matrix); 
for i = 1 : d1 
cycle_index = i; 
Xparam_correction_cycle_row = Xmatrix_obs_rangerate( (cycle_index-1)*Nparam_rangerate_cycle+1 : cycle_index * Nparam_rangerate_cycle, 1  )' ; 
rangerate_emp_param_aposteriori_matrix(cycle_index , Num_time_collumns+1 : d2) = rangerate_emp_param_apriori_matrix(cycle_index , Num_time_collumns+1 : d2) + Xparam_correction_cycle_row ; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Update strucuture array with aposteriori values 
rangerate_empirical_struct.rangerate_empirical_parameters = rangerate_emp_param_aposteriori_matrix; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
end 
grace_pod_struct.rangerate_empirical_struct = rangerate_empirical_struct; 
 
rangerate_postfit_residuals = b_rangerate - A_rangerate * Xmatrix; 
rms_rangerate_postfit_residuals = rms(rangerate_postfit_residuals); 
 
% save rangerate_postfit_residuals.neq  rangerate_postfit_residuals -ASCII -double 
% save rms_rangerate_postfit_residuals.neq  rms_rangerate_postfit_residuals -ASCII -double 
 
