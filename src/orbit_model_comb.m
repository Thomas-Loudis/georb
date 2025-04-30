function [orbit_model_struct] = orbit_model_comb (orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model_comb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Variables for combined parameter estimation 
% Combination of intersatellite ranging observations with GPS-based orbit
% observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               7 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 30/04/2025, Thomas Loudis Papanikolaou
%             Code minor modifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Orbit Configuration structure array
orbit_config_struct = orbit_model_struct.orbit_config;

% Gravity Field structure array
gravity_field_struct = orbit_model_struct.gravity_field;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Parameter estimation settings (TEMP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% Weighted Estimation solution approaches
% 1. Identity matrix 
% 2. Weights based on previous estimation observation residuals
% 3. Fixed weights to sigma defined values
COMBESTIM_weight = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Exclude estimation of gravity field parameters
% 1. Simultaneous parameter estimation of orbit and gravity parameters 
% 2. Step-wise (2 steps) parameter estimation of orbit and gravity parameters 
COMBESTIM_gravity_step = 0;

% Gravity Field structure array
% gravity_field_struct = orbit_model_struct.gravity_field;
gravity_solution_yn = gravity_field_struct.gravity_solution_yn;
test_status_yn = strcmp(gravity_solution_yn,'y');   
if test_status_yn == 1   
COMBESTIM_gravity_step = 1;
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables to structure arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observation_model_matrix.COMBESTIM_combparamestim_01 = COMBESTIM_combparamestim_01;
observation_model_matrix.COMBESTIM_Nestim_comb  = COMBESTIM_Nestim_comb;
observation_model_matrix.COMBESTIM_OBS = COMBESTIM_OBS;

observation_model_matrix.COMBESTIM_intersat_obs = COMBESTIM_intersat_obs;
observation_model_matrix.COMBESTIM_weight = COMBESTIM_weight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gravity_field_determination_matrix.COMBESTIM_gravity_step = COMBESTIM_gravity_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update orbit model structure array
orbit_model_struct.observation_model_matrix             = observation_model_matrix;
orbit_model_struct.gravity_field_determination_matrix   = gravity_field_determination_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
