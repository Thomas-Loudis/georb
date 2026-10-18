function [NEQ_N, NEQ_u, NEQ_N_orbit1, NEQ_u_orbit1, NEQ_N_orbit2, NEQ_u_orbit2, NEQ_N_range, NEQ_u_range, NEQ_N_rangerate, NEQ_u_rangerate] = grace_comb_neqweight(grace_pod_struct, design_matrix_comb, reducedobservation_matrix, NEQ_write, grav_paramestim_yn) 

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
% Function: grace_comb_neqweight 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Normal Equations combination with weights
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
 
 
%--------------------------------------------------------------------------  
% Read data from grace_pod structure array 
orbit_model_GRACE1 = grace_pod_struct.grace1_pod.orbit_model; 
orbit_model_GRACE2 = grace_pod_struct.grace2_pod.orbit_model; 
orbit_model_struct = orbit_model_GRACE1;  
  
% Combined solution Observations type:   
COMBESTIM_OBS = orbit_model_struct.observation_model_matrix.COMBESTIM_OBS;

% Downweight factor of orbits :: applied in combination for gravity solution
downweight_factor_orbit = orbit_model_struct.observation_model_matrix.downweight_factor_orbit;
%--------------------------------------------------------------------------  

%-------------------------------------------------------------------------- 
% Observations Weights 
%-------------------------------------------------------------------------- 
% Weighted estimation solution approaches 
[sigma_obs_orbit1, sigma_obs_orbit2, sigma_intersat_range, sigma_intersat_rangerate] = intersat_weight_function(grace_pod_struct); 
sigma_intersat_ranging = [sigma_intersat_range sigma_intersat_rangerate];

% Gravity Parameter Estimation
test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');    

% DownWeight factor for orbits in case of gravity recovery mode
if test_grav_paramestim_yn == 1    
% downweight_factor_orbit = 1;
sigma_obs_orbit1 = downweight_factor_orbit * sigma_obs_orbit1;
sigma_obs_orbit2 = downweight_factor_orbit * sigma_obs_orbit2;
end
%-------------------------------------------------------------------------- 
% sigma_obs_orbit1
% sigma_obs_orbit2
% sigma_intersat_rangerate

%--------------------------------------------------------------------------  
% Design matrices for Combination solution
A_orbit1    = design_matrix_comb.A_orbit1;
A_orbit2    = design_matrix_comb.A_orbit2;
A_range     = design_matrix_comb.A_range;
A_rangerate = design_matrix_comb.A_rangerate;
Amatrix_param_index = design_matrix_comb.Amatrix_param_index;

%--------------------------------------------------------------------------  
% Reduced-Observation matrices to strucure array
b_orbit1    = reducedobservation_matrix.b_orbit1    ;
b_orbit2    = reducedobservation_matrix.b_orbit2    ;
b_range     = reducedobservation_matrix.b_range     ;
b_rangerate = reducedobservation_matrix.b_rangerate ; 
%--------------------------------------------------------------------------  

%--------------------------------------------------------------------------  
% Normal Equations per block/type of observations 
%--------------------------------------------------------------------------  
% Orbit 1 
% [Xestim_orbit1, NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq_sol(A_orbit1, b_orbit1, sigma_obs_orbit1); 
[NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq(A_orbit1, b_orbit1, sigma_obs_orbit1); 
% Orbit 2 
% [Xestim_orbit2, NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq_sol(A_orbit2, b_orbit2, sigma_obs_orbit2); 
[NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq(A_orbit2, b_orbit2, sigma_obs_orbit2); 
% Range-rate 
% [Xestim_rangerate, NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq_sol(A_rangerate, b_rangerate, sigma_intersat_rangerate); 
[NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq(A_rangerate, b_rangerate, sigma_intersat_rangerate); 
% Range 
% [Xestim_range, NEQ_N_range, NEQ_u_range] = estimator_neq_sol(A_range, b_range, sigma_intersat_range); 
[NEQ_N_range, NEQ_u_range] = estimator_neq(A_range, b_range, sigma_intersat_range); 
%--------------------------------------------------------------------------  

%--------------------------------------------------------------------------  
% Normal Equations combination 
%--------------------------------------------------------------------------  
% Combined solution Observations: Normal Equations sum matrices 
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
%--------------------------------------------------------------------------  
 

