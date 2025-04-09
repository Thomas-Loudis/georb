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


orbit_config_struct = orbit_model_struct.orbit_config_struct;

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

% Combined solution approach to Gravity and Orbit parameters estimation:  
% 1. Combined estimation : Gravity and Orbit parameters combined estimation
% 4. Gravity parameters estimation only - Orbit parameters are fixed to apriori POD solution
%COMBESTIM_comb_solution = 4
COMBESTIM_comb_solution = 1;

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
COMBESTIM_intersat_obs_LRI = 1;

% Weighted Estimation solution approaches
% 1. Identity matrix 
% 2. Weights based on previous estimation observation residuals
% 3. Fixed weights to sigma defined values
COMBESTIM_weight = 3;

% Integration Step size reduced for LRI/KBR observations combined estimation
COMBESTIM_integration_step = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Exclude estimation of gravity field parameters
% 1. Simultaneous parameter estimation of orbit and gravity parameters 
% 2. Step-wise (2 steps) parameter estimation of orbit and gravity parameters 
COMBESTIM_gravity_step = 0;

% Gravity parmeters partials (in case COMBESTIM_gravity_step = 0)
COMBESTIM_gravity_partials = 1;
COMBESTIM_gravity_partials = 0;

% Re-read orbit_model when update configuration file for gravity parameters
paramestim_update_config_01 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables to structure arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observation_model_matrix.COMBESTIM_combparamestim_01 = COMBESTIM_combparamestim_01;
observation_model_matrix.COMBESTIM_Nestim_comb  = COMBESTIM_Nestim_comb;
observation_model_matrix.COMBESTIM_comb_solution = COMBESTIM_comb_solution;
observation_model_matrix.COMBESTIM_OBS = COMBESTIM_OBS;

observation_model_matrix.COMBESTIM_intersat_obs_LRI = COMBESTIM_intersat_obs_LRI;
observation_model_matrix.COMBESTIM_weight = COMBESTIM_weight;
observation_model_matrix.COMBESTIM_integration_step = COMBESTIM_integration_step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gravity_field_determination_matrix.COMBESTIM_gravity_step = COMBESTIM_gravity_step;
gravity_field_determination_matrix.COMBESTIM_gravity_partials = COMBESTIM_gravity_partials;

gravity_field_determination_matrix.paramestim_update_config_01 = paramestim_update_config_01;

% gravity_model_matrix = gravity_field_determination_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update orbit model structure array
orbit_model_struct.observation_model_matrix = observation_model_matrix;

orbit_model_struct.gravity_field_determination_matrix = gravity_field_determination_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
