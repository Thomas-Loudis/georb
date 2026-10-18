function [orbit_model_struct] = orbit_model_comb (orbit_model_struct) 

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
% Function: orbit_model_comb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Set Variables for combined parameter estimation  
%  Combination of intersatellite ranging observations with GPS-based orbit 
%  observations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                               7 April 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 30/04/2025, Thomas Loudis Papanikolaou 
%             Code minor modifications  
% 10/06/2026, Thomas Loudis Papanikolaou 
%             Code minor modifications  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
% Orbit Configuration structure array 
orbit_config_struct = orbit_model_struct.orbit_config; 
   
%-------------------------------------------------------------------------- 
% Combined Parameter estimation settings 
%-------------------------------------------------------------------------- 
% y/n status 
COMBESTIM_combparamestim_01 = 0; 
param_keyword = 'combined_param_estim_yn'; 
[combined_param_estim_yn] = read_param_cfg(orbit_config_struct,param_keyword); 
test_status_yn = strcmp(combined_param_estim_yn,'y');    
if test_status_yn == 1  
COMBESTIM_combparamestim_01 = 1; 
end 
 
% Number of estimation occurences (If 0 :: Combined-Estimator not applied)    
COMBESTIM_Nestim_comb = 1; 
 
% Combined solution Observations:   
% 1. OBS: Orbit1 
% 2. OBS: Orbit2 
% 3. OBS: Orbit1 + Orbit2 
% 4. OBS: LRI range 
% 5. OBS: LRI rangerate 
% 6. OBS: LRI range + rangerate 
% 7. OBS: Orbit1 + Orbit2 + LRI range 
% 8. OBS: Orbit1 + Orbit2 + LRI range + rangerate 
% 9. OBS: Orbit1 + Orbit2 + LRI rangerate 
COMBESTIM_OBS = 9; 
 
% Intersatellite Observations :: Laser (LRI) or K-band (KBR) 
% 0. KBR: COMBESTIM_intersat_obs_LRI = 0 
% 1. LRI: COMBESTIM_intersat_obs_LRI = 1 
% COMBESTIM_intersat_obs_LRI = 1; 
param_keyword = 'KBR_obs_estim'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
test_status_yn = strcmp(param_value,'y');    
if test_status_yn == 1  
COMBESTIM_intersat_obs = 'KBR'; 
end 
 
param_keyword = 'LRI_obs_estim'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
test_status_yn = strcmp(param_value,'y');    
if test_status_yn == 1  
COMBESTIM_intersat_obs = 'LRI'; 
end 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Weighted Estimation solution approaches 
%-------------------------------------------------------------------------- 
% Weights approach of combined parameter estimator
% 1. Identity matrix  
% 2. VCE: Weights based on observation covariance matrix Cy of previous estimation  
% 3. Fixed weights to sigma defined values 
% 4. Weights based on observation residuals of previous estimation 
param_keyword = 'combined_weights_variances_yn'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
combined_weights_variances_yn = param_value;
test_status_yn = strcmp(combined_weights_variances_yn,'y');    
if test_status_yn == 1  
    COMBESTIM_weight = 2;
else
% Fixed Weights :: Case 3
    COMBESTIM_weight = 3;
end 

% Fixed Weights case 
% Weights are computed based on the a-priori Sigma values of the orbits observations, range and range-rate observations 
% Sigma Units are cm, mm, μm/sec for orbits, inter-satellite range, inter-satellite range-rate respectively   
param_keyword = 'sigma_obs_orbit'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
sigma_obs_orbit = sscanf(param_value,'%e %*'); 
 
param_keyword = 'sigma_intersat_range'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
sigma_intersat_range = sscanf(param_value,'%e %*'); 
 
param_keyword = 'sigma_intersat_rangerate'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
sigma_intersat_rangerate = sscanf(param_value,'%e %*'); 

param_keyword = 'downweight_factor_orbit'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
downweight_factor_orbit = sscanf(param_value,'%d %*'); 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Case for GRACE-like missions (GRACE, MAGIC, NGGM) with accelerometers
% Accelerometers Scale factors consideration as global parameters or not (y/n) ::
% - y: Global/common parameters along with gravity parameters. They will be included in the reduced Normal Equations matrices data (arc/daily NEQs) following the pre-elimination of the orbital arc-based parameters) 
% - n: Orbit arc parameters (excluded from the reduced Normal Equations matrices)
% acc_scale_global_yn           y

param_keyword = 'acc_scale_global_yn'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
acc_scale_global_yn = param_value;

% Update Gravity Field structure array 
orbit_model_struct.gravity_field.acc_scale_global_yn = acc_scale_global_yn; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Intersatellite Range-Rate empirical parameters
%-------------------------------------------------------------------------- 
% y/n status 
param_keyword = 'intersat_rangerate_empirical_yn'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
intersat_rangerate_empirical_yn = param_value;

param_keyword = 'intersat_rangerate_empirical_step'; 
[param_value] = read_param_cfg(orbit_config_struct,param_keyword); 
intersat_rangerate_empirical_step = sscanf(param_value,'%e %*');

% structure array
rangerate_empirical_struct.intersat_rangerate_empirical_yn   = intersat_rangerate_empirical_yn;
rangerate_empirical_struct.intersat_rangerate_empirical_step = intersat_rangerate_empirical_step;
%-------------------------------------------------------------------------- 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Gravity Field parameter estimation approach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Gravity Field structure array 
% gravity_field_struct = orbit_model_struct.gravity_field; 
% gravity_solution_yn = gravity_field_struct.gravity_solution_yn; 
% test_status_yn = strcmp(gravity_solution_yn,'y');    
% COMBESTIM_gravity_step = 0; 
% if test_status_yn == 1    
% COMBESTIM_gravity_step = 1; 
% end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Variables to structure arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
observation_model_matrix.COMBESTIM_combparamestim_01    = COMBESTIM_combparamestim_01; 
observation_model_matrix.COMBESTIM_Nestim_comb          = COMBESTIM_Nestim_comb; 
observation_model_matrix.COMBESTIM_OBS                  = COMBESTIM_OBS; 
 
observation_model_matrix.COMBESTIM_intersat_obs = COMBESTIM_intersat_obs; 
observation_model_matrix.COMBESTIM_weight       = COMBESTIM_weight; 
 
observation_model_matrix.sigma_obs_orbit          = sigma_obs_orbit; 
observation_model_matrix.sigma_intersat_range     = sigma_intersat_range; 
observation_model_matrix.sigma_intersat_rangerate = sigma_intersat_rangerate; 

observation_model_matrix.downweight_factor_orbit = downweight_factor_orbit; 

% gravity_field_determination_matrix.COMBESTIM_gravity_step = COMBESTIM_gravity_step; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Update orbit model structure array 
orbit_model_struct.observation_model_matrix             = observation_model_matrix; 
% orbit_model_struct.gravity_field_determination_matrix   = gravity_field_determination_matrix; 
orbit_model_struct.rangerate_empirical_struct           = rangerate_empirical_struct; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
