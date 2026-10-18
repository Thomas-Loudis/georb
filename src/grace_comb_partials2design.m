function [design_matrix_perobs, reducedobservation_matrix, Nparam_orbital] = grace_comb_partials2design(grace_pod_struct, Nglobal_grav, Nglobal_scale, rangerate_empirical_struct) 

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
% Function: grace_comb_partials2design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Processing of Partials matrices from Variational Equations solution to
%  form the design matrices per observation type (orbits and
%  inter-satellite ranging observations).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbit_model_struct  : Orbit model structure array  
% - ic_data_struct      : Initial Conditions structure array  
% 
% Output arguments: 
% - out_dir_name        : Output Data folder name  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                          28 September 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Read data from grace_pod structure array 
orbit_model_GRACE1 = grace_pod_struct.grace1_pod.orbit_model; 
orbit_model_GRACE2 = grace_pod_struct.grace2_pod.orbit_model; 
orbit_model_struct = orbit_model_GRACE1;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
% Combined solution Observations type:   
COMBESTIM_OBS = orbit_model_struct.observation_model_matrix.COMBESTIM_OBS; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Read data from grace_pod structure array 
orbcGA                       = grace_pod_struct.grace1_pod.orbit_matrices.orbit_crf; 
orbcGB                       = grace_pod_struct.grace2_pod.orbit_matrices.orbit_crf; 
sat1_OBS_matrix              = grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix; 
sat2_OBS_matrix              = grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix; 
 
% sat1_OBS_residuals           = grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals; 
% sat2_OBS_residuals           = grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals; 
 
sat1_veqZ_matrix             = grace_pod_struct.grace1_pod.orbit_matrices.state_transition_matrix; 
sat1_veqP_matrix             = grace_pod_struct.grace1_pod.orbit_matrices.sensitivity_matrix; 
 
sat2_veqZ_matrix             = grace_pod_struct.grace2_pod.orbit_matrices.state_transition_matrix; 
sat2_veqP_matrix             = grace_pod_struct.grace2_pod.orbit_matrices.sensitivity_matrix; 
 
intersat_observation_data    = grace_pod_struct.intersat_pod.intersat_observation_data; 
% intersat_range_residuals     = grace_pod_struct.intersat_pod.intersat_range_residuals; 
% intersat_rangerate_residuals = grace_pod_struct.intersat_pod.intersat_rangerate_residuals; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%-------------------------------------------------------------------------- 
% Range-Rate Empirical parameters 
%-------------------------------------------------------------------------- 
rangerate_emp_01       = rangerate_empirical_struct.rangerate_empirical_function_01;
rangerate_emp_param    = rangerate_empirical_struct.rangerate_empirical_parameters; 
Nparam_rangerate_cycle = rangerate_empirical_struct.Nparam_rangerate_cycle;   
Nparam_rangerate_all   = rangerate_empirical_struct.Nparam_rangerate ;  
if rangerate_emp_01 == 1 
Nparam_obs_rangerate = Nparam_rangerate_all; 
else 
Nparam_obs_rangerate = 0;     
end 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% sigma prior combination (NO weights are needed at this step)
sigma_obs_orbit1 = 1;
sigma_obs_orbit2 = 1;
sigma_intersat_rangerate = 1;
sigma_intersat_range = 1;
sigma_intersat_ranging = [sigma_intersat_range sigma_intersat_rangerate];
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% GRACE-4 Design and NEQ matrices  
%-------------------------------------------------------------------------- 
% Inter-satellite observations w.r.t. GRACE-2/-4 
[Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = neq_intersat_ranging(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, sigma_intersat_ranging, rangerate_empirical_struct);   
% Orbit pseudo-observations 
[NEQn_grace2, NEQu_grace2, Amatrix_obsorb_grace2, Wmatrix_obsorb_grace2] = neq_orbit (orbcGB, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, sigma_obs_orbit2); 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% GRACE-3 Design and NEQ matrices  
%-------------------------------------------------------------------------- 
% Inter-satellite observations w.r.t. GRACE-1/-3 
[Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = neq_intersat_ranging(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, sigma_intersat_ranging, rangerate_empirical_struct);  
% Orbit pseudo-observations 
[NEQn_grace1, NEQu_grace1, Amatrix_obsorb_grace1, Wmatrix_obsorb_grace1] = neq_orbit (orbcGA, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, sigma_obs_orbit1); 
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------   
% Design matrices per observation type
design_matrix_perobs.rangerate_sat1    = Amatrix_rangerate_sat1;
design_matrix_perobs.range_sat1        = Amatrix_range_sat1;
design_matrix_perobs.rangerate_sat2    = Amatrix_rangerate_sat2;
design_matrix_perobs.range_sat2        = Amatrix_range_sat2;
design_matrix_perobs.obsorb_grace1     = Amatrix_obsorb_grace1;
design_matrix_perobs.obsorb_grace2     = Amatrix_obsorb_grace2;
%--------------------------------------------------------------------------   

%--------------------------------------------------------------------------   
% Reduced-Observations matrix 
%--------------------------------------------------------------------------   
b_orbit1    = Wmatrix_obsorb_grace1; 
b_orbit2    = Wmatrix_obsorb_grace2; 
b_range     = Wmatrix_range_sat2; 
b_rangerate = Wmatrix_rangerate_sat2; 

% Reduced-Observation matrices to strucure array
reducedobservation_matrix.b_orbit1    = b_orbit1;
reducedobservation_matrix.b_orbit2    = b_orbit2;
reducedobservation_matrix.b_range     = b_range;
reducedobservation_matrix.b_rangerate = b_rangerate;
%--------------------------------------------------------------------------   

%--------------------------------------------------------------------------   
% Number of parameters, observations  
%--------------------------------------------------------------------------   
[d1, d2] = size(Amatrix_obsorb_grace1); 
[d3, d4] = size(Amatrix_obsorb_grace2); 
Nparam_grace1 = d2; 
Nparam_grace2 = d4; 
Nobs_sat1 = d1; 
Nobs_sat2 = d3; 
% Overall parameters number of the combined parameter estimation including
% gravity parameters (global parameters) 
Nparam_COMB = (Nparam_grace1 - Nglobal_grav) + (Nparam_grace2 - Nglobal_grav) + Nglobal_grav + Nparam_obs_rangerate;
% Number of orbital parameters (individual arc parameters) 
Nparam_orbital = Nparam_COMB - Nglobal_grav;
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------   
% Design Matrix combination :: Design matrix level processing for forming the NEQ combination  
%--------------------------------------------------------------------------   
% % Design matrix combination
% [A_orbit1, A_orbit2, A_range, A_rangerate, Amatrix_param_index] = grace_comb_design(design_matrix_perobs, Nglobal_grav, Nglobal_scale, Nparam_obs_rangerate);
% 
% % Design matrices for Combination solution
% design_matrix_comb.A_orbit1    = A_orbit1;
% design_matrix_comb.A_orbit2    = A_orbit2;
% design_matrix_comb.A_range     = A_range;
% design_matrix_comb.A_rangerate = A_rangerate;
% design_matrix_comb.Amatrix_param_index = Amatrix_param_index;
%--------------------------------------------------------------------------   
