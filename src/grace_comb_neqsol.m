function [Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xgrav, grace_pod_struct] = grace_comb_neqsol(grace_pod_struct, NEQ_N, NEQ_u, design_matrix_comb, reducedobservation_matrix, rangerate_empirical_struct, gravity_arc_sol_flag, Nglobal_scale) 

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
% 28/09/2026, Thomas Papanikolaou 
%             Source Code minor modification   
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Range-Rate Empirical parameters 
%-------------------------------------------------------------------------- 
% Parameter estimation On/Off 1/0
rangerate_emp_01       = rangerate_empirical_struct.rangerate_empirical_function_01;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Orbit model
orbit_model_G1  = grace_pod_struct.grace1_pod.orbit_model; 
orbit_model_G2  = grace_pod_struct.grace2_pod.orbit_model; 
orbit_model_struct = orbit_model_G1;
% Combined solution Observations type:   
COMBESTIM_OBS = orbit_model_struct.observation_model_matrix.COMBESTIM_OBS; 
%-------------------------------------------------------------------------- 


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
% Parameter partials collumns indexes
%-------------------------------------------------------------------------- 
col_orbit1_1 = Amatrix_param_index.orbit1(1,1);
col_orbit1_2 = Amatrix_param_index.orbit1(1,2);

col_orbit2_1 = Amatrix_param_index.orbit2(1,1);
col_orbit2_2 = Amatrix_param_index.orbit2(1,2);

% Gravity parameters
col_common_1 = Amatrix_param_index.global(1,1);
col_common_2 = Amatrix_param_index.global(1,2);

col_obs_rangerate_1 = Amatrix_param_index.rangerate_emp(1,1);
col_obs_rangerate_2 = Amatrix_param_index.rangerate_emp(1,2);

col_orbit1_ACC_Scale_start  = Amatrix_param_index.grace1_acc_scale(1,1);
col_orbit1_ACC_Scale_end    = Amatrix_param_index.grace1_acc_scale(1,2);

col_orbit2_ACC_Scale_start  = Amatrix_param_index.grace2_acc_scale(1,1);
col_orbit2_ACC_Scale_end    = Amatrix_param_index.grace2_acc_scale(1,2);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Combined solution Observations: 
if COMBESTIM_OBS == 1 
b_comb = b_orbit1; 
A_comb = A_orbit1; 
 
elseif COMBESTIM_OBS == 2 
b_comb = b_orbit2; 
A_comb = A_orbit2; 
     
elseif COMBESTIM_OBS == 3 
b_comb = [b_orbit1; b_orbit2]; 
A_comb = [A_orbit1; A_orbit2];
     
elseif COMBESTIM_OBS == 4 
b_comb = [b_range]; 
A_comb = [A_range];
   
elseif COMBESTIM_OBS == 5 
b_comb = [b_rangerate]; 
A_comb = [A_rangerate];
 
elseif COMBESTIM_OBS == 6 
b_comb = [b_range; b_rangerate]; 
A_comb = [A_range; A_rangerate];
 
elseif COMBESTIM_OBS == 7  
b_comb = [b_orbit1; b_orbit2; b_range]; 
A_comb = [A_orbit1; A_orbit2; A_range];
 
elseif COMBESTIM_OBS == 8 
b_comb = [b_orbit1; b_orbit2; b_range; b_rangerate]; 
A_comb = [A_orbit1; A_orbit2; A_range; A_rangerate];

elseif COMBESTIM_OBS == 9 
b_comb = [b_orbit1; b_orbit2; b_rangerate]; 
A_comb = [A_orbit1; A_orbit2; A_rangerate];
end 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Normal Equations Solution
%-------------------------------------------------------------------------- 
% Matrix inversion algorithm ID no. :: Cholesky factorisation 
inv_id = 6;

% NEQ scaling tolerance    
NEQscaling_condtol = 10^7;

% Normal Equations solution
[Xmatrix, NEQ_N_scal, NEQ_u_scal, NEQscaling_mat] = neq_sol(NEQ_N, NEQ_u, inv_id, NEQscaling_condtol);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Reference Variance
%-------------------------------------------------------------------------- 
m_param = size(Xmatrix,1) ;
n_obs = size(b_comb,1) ;
error_matrix = b_comb - A_comb * Xmatrix;
error_matrix_square = (error_matrix' * error_matrix) ;
sigma0 = sqrt( error_matrix_square / (n_obs - m_param) ); 
% sigma0_grace_comb_neqsol = sigma0
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Covariance matrices 
%-------------------------------------------------------------------------- 
% Covariance matrix of parameters 
[Nmatrix_inv] = inv_mat(NEQ_N_scal, inv_id); 
Cx_matrix = sigma0^2 * Nmatrix_inv;
% Estimated Parameter de-scaling (in case of non-scaling NEQscaling_mat=1)
Cx_matrix = NEQscaling_mat * Cx_matrix * NEQscaling_mat; 
sigma_Xparam_combined = sqrt( diag(Cx_matrix) ); 
% save sigma_Xparam_combined.neq  sigma_Xparam_combined -ASCII -double 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Orbit Parameters estimated corrections     
%--------------------------------------------------------------------------
if Nglobal_scale == 0
Xmatrix_orbit1 = Xmatrix(1 : col_orbit1_2 , 1); 
Xmatrix_orbit2 = Xmatrix(col_orbit2_1 : col_orbit2_2 , 1) ;  
elseif Nglobal_scale > 0
Xmatrix_orbit1 = [Xmatrix(1 : col_orbit1_2 , 1); Xmatrix(col_orbit1_ACC_Scale_start : col_orbit1_ACC_Scale_end , 1)] ; 
Xmatrix_orbit2 = [Xmatrix(col_orbit2_1 : col_orbit2_2 , 1); Xmatrix(col_orbit2_ACC_Scale_start : col_orbit2_ACC_Scale_end , 1)] ;  
end

% Structure array update :: param_aposteriori 
grace_pod_struct.grace1_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit1; 
grace_pod_struct.grace2_pod.orbit_matrices.param_cor_Xmatrix = Xmatrix_orbit2; 
 
grace_pod_struct.grace1_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit1; 
grace_pod_struct.grace2_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit2; 
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Gravity Field parameters estimated corrections 
%-------------------------------------------------------------------------- 
if gravity_arc_sol_flag == 1 
Xmatrix_common = Xmatrix(col_common_1 : col_common_2 , 1); 
Xcommon        = Xmatrix_common; 
else 
Xcommon = 0; 
end
Xgrav = Xcommon;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Range-Rate Empirical Parameters matrix :: Add estimated corrections 
%-------------------------------------------------------------------------- 
if rangerate_emp_01 == 1 
Xmatrix_obs_rangerate = Xmatrix(col_obs_rangerate_1 : col_obs_rangerate_2 , 1); 
   
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
 
% Strucuture array: assign aposteriori values 
rangerate_empirical_struct.rangerate_empirical_parameters = rangerate_emp_param_aposteriori_matrix; 

end
% Structure array update :: Range-Rate Empirical Parameters  
grace_pod_struct.rangerate_empirical_struct = rangerate_empirical_struct; 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Residuals Statistics
%-------------------------------------------------------------------------- 
% Range-Rate
rangerate_postfit_residuals = b_rangerate - A_rangerate * Xmatrix; 
rms_rangerate_postfit_residuals = rms(rangerate_postfit_residuals);

% Range
range_postfit_residuals = b_range - A_range * Xmatrix; 
rms_range_postfit_residuals = rms(range_postfit_residuals); 

% Orbit-1
orbit1_postfit_residuals = b_orbit1 - A_orbit1 * Xmatrix; 
% orbit1_postfit_residuals = b_orbit1 - Amatrix_obsorb_grace1 * Xmatrix_orbit1_scaled; 
rms_orbit1_postfit_residuals = rms(orbit1_postfit_residuals);

% Orbit-2
orbit2_postfit_residuals = b_orbit2 - A_orbit2 * Xmatrix; 
rms_orbit2_postfit_residuals = rms(orbit2_postfit_residuals);
%-------------------------------------------------------------------------- 
% save rangerate_postfit_residuals.neq  rangerate_postfit_residuals -ASCII -double 
% save rms_rangerate_postfit_residuals.neq  rms_rangerate_postfit_residuals -ASCII -double 


%-------------------------------------------------------------------------- 
% Structure array update :: InterSatellite Range and Range-Rate Residuals  
%-------------------------------------------------------------------------- 
% Number of time arguments 
Ntime_col = 1; 

% Apriori Inter-Satellite ranging residuals
intersat_range_residuals_apriori = grace_pod_struct.intersat_pod.intersat_range_residuals;
intersat_rangerate_residuals_apriori = grace_pod_struct.intersat_pod.intersat_rangerate_residuals;
% Epochs
intersat_range_epochs       = intersat_range_residuals_apriori(:,1);
intersat_rangerate_epochs   = intersat_rangerate_residuals_apriori(:,1);

grace_pod_struct.intersat_pod.intersat_range_residuals      = [intersat_range_epochs range_postfit_residuals];
grace_pod_struct.intersat_pod.intersat_rangerate_residuals  = [intersat_rangerate_epochs rangerate_postfit_residuals];
%-------------------------------------------------------------------------- 

if 1>0 
fprintf('%s \n', 'Post-Fit residuals:'); 
fprintf('%s \n', 'Intersatellite Ranging Post-Fit residuals:'); 
fprintf('%s%17.6f', 'Range-Rate residuals: RMS (μm/sec): ',rms_rangerate_postfit_residuals * 10^6); 
fprintf('\n'); 
fprintf('%s%17.6f', 'Range residuals: RMS (mm):          ',rms_range_postfit_residuals * 10^3); 
fprintf('\n'); 

fprintf('%s%17.6f', 'Orbit1 residuals: RMS (m):          ',rms_orbit1_postfit_residuals); 
fprintf('\n'); 
fprintf('%s%17.6f', 'Orbit2 residuals: RMS (m):          ',rms_orbit2_postfit_residuals); 
fprintf('\n'); 
fprintf('\n') 
end 

%-------------------------------------------------------------------------- 
% Reference Variances per individual observation set
%-------------------------------------------------------------------------- 
m_param = size(Xmatrix,1);

% Range-Rate
n_obs = size(b_rangerate,1);
error_matrix = rangerate_postfit_residuals;
sigma = sqrt( (error_matrix' * error_matrix) / (n_obs - m_param) ); 
sigma0_intersat_rangerate = sigma ; 

% Rang
n_obs = size(b_range,1);
error_matrix = range_postfit_residuals;
sigma = sqrt( (error_matrix' * error_matrix) / (n_obs - m_param) ); 
sigma0_intersat_range = sigma ;

% Orbit-1
n_obs = size(b_orbit1,1);
error_matrix = orbit1_postfit_residuals;
sigma = sqrt( (error_matrix' * error_matrix) / (n_obs - m_param) ); 
sigma0_orbit1 = sigma ;
Cy_orbit1 = sigma^2;

% Orbit-2
n_obs = size(b_orbit2,1);
error_matrix = orbit2_postfit_residuals;
sigma = sqrt( (error_matrix' * error_matrix) / (n_obs - m_param) ); 
sigma0_orbit2 = sigma ;
Cy_orbit2 = sigma^2;
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% grace_pod_struct update :: overall structure array :: Assign updated values to matrices 
%-------------------------------------------------------------------------- 
% Residuals
grace_pod_struct.intersat_pod.intersat_range_residuals      = range_postfit_residuals; 
grace_pod_struct.intersat_pod.intersat_rangerate_residuals  = rangerate_postfit_residuals; 

% grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals = orbit1_postfit_residuals; 
% grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals = orbit2_postfit_residuals; 

% grace_pod_struct.grace1_pod.orbit_matrices.rms_observation_residuals = rms_orbit1_postfit_residuals; 
% grace_pod_struct.grace2_pod.orbit_matrices.rms_observation_residuals = rms_orbit2_postfit_residuals; 


% Statistics Variances Components
grace_pod_struct.grace1_pod.orbit_model.sigma0 = sigma0_orbit1; 
grace_pod_struct.grace2_pod.orbit_model.sigma0 = sigma0_orbit2;

grace_pod_struct.grace1_pod.orbit_model.Cy = Cy_orbit1; 
grace_pod_struct.grace2_pod.orbit_model.Cy = Cy_orbit2; 
%-------------------------------------------------------------------------- 
