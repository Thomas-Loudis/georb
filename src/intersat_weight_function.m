function [sigma_obs_orbit1, sigma_obs_orbit2, sigma_intersat_range, sigma_intersat_rangerate] = intersat_weight_function(grace_pod_struct) 

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
% - grace_pod_struct            : Overall structure array output of grace_pod function  
% 
% Output arguments: 
% - sigma_obs_orbit1            : Weight for orbit observations of GRACE satellite 1  
% - sigma_obs_orbit2            : Weight for orbit observations of GRACE satellite 2  
% - sigma_intersat_range        : Weight for inter-satellite range observations  
% - sigma_intersat_rangerate    : Weight for inter-satellite range-rate observations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                              22 April 2025 
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
% Combined Parameter estimation settings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
observation_model_matrix = orbit_model_struct.observation_model_matrix; 
 
% Weighted Estimation solution approaches 
COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
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
% Observation Weights 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Weighted estimation solution approaches 
weight_sol_opt = COMBESTIM_weight; 
     
% Weights based on errors/residuals from step-1 orbit parameter estimation   
if weight_sol_opt == 4 
% Orbit Residuals matrices 
[n1, n2] = size(sat1_OBS_residuals); 
for i = 1 : n1 
sat1_OBS_residuals_vector( (i-1)*3+1 : 3*i , 1) = [sat1_OBS_residuals(i,2); sat1_OBS_residuals(i,3); sat1_OBS_residuals(i,4)]'; 
end 
 
[n1, n2] = size(sat2_OBS_residuals); 
for i = 1 : n1 
sat2_OBS_residuals_vector( (i-1)*3+1 : 3*i , 1) = [sat2_OBS_residuals(i,2); sat2_OBS_residuals(i,3); sat2_OBS_residuals(i,4)]'; 
end 
 
[RMS_orbit1_residuals, sigma_orbit1_residuals] = rms(sat1_OBS_residuals_vector); 
[RMS_orbit2_residuals, sigma_orbit2_residuals] = rms(sat2_OBS_residuals_vector); 
 
% Intersatellite (LRI) range and range-rate residuals 
[RMS_range_residuals, sigma_range_residuals] = rms(intersat_range_residuals(:,2)); 
[RMS_rangerate_residuals, sigma_rangerate_residuals] = rms(intersat_rangerate_residuals(:,2)); 
 
% Output 
sigma_obs_orbit1 = sigma_orbit1_residuals; 
sigma_obs_orbit2 = sigma_orbit2_residuals; 
sigma_intersat_range = sigma_range_residuals; 
sigma_intersat_rangerate = sigma_rangerate_residuals; 
 
 
% Weights based on covariance matrices from step-1 orbit parameter estimation   
elseif weight_sol_opt == 2 
Cv_sat1_est = orbit_model_GRACE1.Cy; 
Cv_sat2_est = orbit_model_GRACE2.Cy; 
 
% Sigma factor to the orbits' Cv_sqrt 
% sigma_obsorb = 2 * 10^-2; 
% sigma_obsorb = 10  
sigma_obsorb = 5; 
% sigma_obsorb = 100  
 
% Output 
sigma_obs_orbit1 = sigma_obsorb * sqrt(Cv_sat1_est); 
sigma_obs_orbit2 = sigma_obsorb * sqrt(Cv_sat2_est); 
 
sigma_range = 2 * 10^-3; 
sigma_rangerate = 0.2 * 10^-6; 
 
% Output 
sigma_intersat_range = sigma_range; 
sigma_intersat_rangerate = sigma_rangerate; 
 
 
% Fixed weights to sigma defined values 
elseif weight_sol_opt == 3 
sigma_obs_orbit_config = observation_model_matrix.sigma_obs_orbit; 
sigma_intersat_range_config = observation_model_matrix.sigma_intersat_range; 
sigma_intersat_rangerate_config = observation_model_matrix.sigma_intersat_rangerate; 
 
sigma_obsorb    = sigma_obs_orbit_config * 10^-2; 
sigma_range     = sigma_intersat_range_config * 10^-3; 
sigma_rangerate = sigma_intersat_rangerate_config * 10^-6; 
 
% Output 
sigma_obs_orbit1 = sigma_obsorb; 
sigma_obs_orbit2 = sigma_obsorb; 
sigma_intersat_range = sigma_range; 
sigma_intersat_rangerate = sigma_rangerate; 
 
 
% Identity matrix 
elseif weight_sol_opt == 1 
sigma_obs_orbit1 = 1; 
sigma_obs_orbit2 = 1; 
sigma_intersat_range = 1; 
sigma_intersat_rangerate = 1; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% save sigma_obs_orbit2.neq  sigma_obs_orbit2 -ASCII -double 
% save sigma_obs_orbit1.neq  sigma_obs_orbit1 -ASCII -double 
