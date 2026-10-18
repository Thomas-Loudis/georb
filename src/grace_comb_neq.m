function [NEQ_N, NEQ_u, NEQ_N_reduced, NEQ_u_reduced, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xgrav, sigma_Xgrav, grace_pod_struct] = grace_comb_neq(grace_pod_struct, acc_scale_common_01, NEQ_write, gravity_arc_sol_flag) 

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
%             Source Code upgrade and introduction of new functions   
%-------------------------------------------------------------------------- 
 

%-------------------------------------------------------------------------- 
% Global Parameters Number
% Gravity field coefficients and (optional) accelerometer scale factors
%-------------------------------------------------------------------------- 
[Nglobal_grav, Nglobal_scale, grav_paramestim_yn] = grace_comb_globalparam(grace_pod_struct, acc_scale_common_01);
N_param_GRAV = Nglobal_grav;
N_common_ACC_Scale = Nglobal_scale;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Range-Rate Empirical parameters 
%-------------------------------------------------------------------------- 
% Range-Rate Empirical Parameters :: Structure array 
rangerate_empirical_struct = grace_pod_struct.rangerate_empirical_struct ; 

rangerate_emp_01       = rangerate_empirical_struct.rangerate_empirical_function_01;
Nparam_rangerate_all   = rangerate_empirical_struct.Nparam_rangerate ;  
% Number of Range-Rate empirial parameters
Nparam_obs_rangerate = Nparam_rangerate_all * rangerate_emp_01; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Partials to Design Matrices and Reduced-Observation matrices for combination
% [design_matrix_comb, design_matrix_perobs, reducedobservation_matrix, Nparam_orbital] = grace_comb_partials2design(grace_pod_struct, Nglobal_grav, Nglobal_scale, rangerate_empirical_struct) ;
[design_matrix_perobs, reducedobservation_matrix, Nparam_orbital] = grace_comb_partials2design(grace_pod_struct, Nglobal_grav, Nglobal_scale, rangerate_empirical_struct) ;
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------   
% Design Matrix combination :: Design matrix level processing for forming the NEQ combination  
%--------------------------------------------------------------------------   
% Design matrix combination
[A_orbit1, A_orbit2, A_range, A_rangerate, Amatrix_param_index] = grace_comb_design(design_matrix_perobs, Nglobal_grav, Nglobal_scale, Nparam_obs_rangerate);
% Design matrices for Combination solution
design_matrix_comb.A_orbit1    = A_orbit1;
design_matrix_comb.A_orbit2    = A_orbit2;
design_matrix_comb.A_range     = A_range;
design_matrix_comb.A_rangerate = A_rangerate;
design_matrix_comb.Amatrix_param_index = Amatrix_param_index;
%--------------------------------------------------------------------------   

%--------------------------------------------------------------------------  
% Normal Equations combination 
%--------------------------------------------------------------------------  
[NEQ_N, NEQ_u, NEQ_N_orbit1, NEQ_u_orbit1, NEQ_N_orbit2, NEQ_u_orbit2, NEQ_N_range, NEQ_u_range, NEQ_N_rangerate, NEQ_u_rangerate] = grace_comb_neqweight(grace_pod_struct, design_matrix_comb, reducedobservation_matrix, NEQ_write, grav_paramestim_yn);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Gravity Parameter estimation test
test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Reduced NEQ matrices :: pre-elimination of orbital arc-related parameters 
%-------------------------------------------------------------------------- 
if test_grav_paramestim_yn == 1            
% n_orbparam  = (Nparam_grace1 - N_param_GRAV) + (Nparam_grace2 - N_param_GRAV) + Nparam_obs_rangerate - N_common_ACC_Scale; 
n_orbparam    = Nparam_orbital - N_common_ACC_Scale;
n_globalparam = N_param_GRAV   + N_common_ACC_Scale;
[NEQ_N_reduced, NEQ_u_reduced, Nzz, Nzg, Ngg, u_z, u_g] = neq_reduced(NEQ_N, NEQ_u, n_orbparam, n_globalparam);
else 
NEQ_N_reduced = 0; 
NEQ_u_reduced = 0; 
end 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Gravity Field solution based on orbit arc NEQs following orbit parameters pre-elimination  
%-------------------------------------------------------------------------- 
if test_grav_paramestim_yn == 1 && gravity_arc_sol_flag == 1
% Matrix inversion algorithm ID no.
inv_id = 6;
% NEQ scaling tolerance    
NEQscaling_condtol = 10^7;
% Normal Equations solution
[Xmatrix_NEQ_reduced, NEQ_N_reduced_scal, NEQ_u_reduced_scal, NEQscaling_mat] = neq_sol(NEQ_N_reduced, NEQ_u_reduced, inv_id, NEQscaling_condtol);
Xgrav_neqreduced = Xmatrix_NEQ_reduced(N_common_ACC_Scale +1 : end,1);

%-------------------------------------------------------------------------- 
% Covariance matrix 
%--------------------------------------------------------------------------
NEQ_N_scal = NEQ_N_reduced_scal;
sigma0 = 1;  
% Covariance matrix of parameters 
[Nmatrix_inv] = inv_mat(NEQ_N_scal, inv_id); 
Cx_matrix = sigma0^2 * Nmatrix_inv;
% Estimated Parameter de-scaling (in case of non-scaling NEQscaling_mat=1)
Cx_matrix = NEQscaling_mat * Cx_matrix * NEQscaling_mat; 
sigma_x = sqrt( diag(Cx_matrix) ); 
sigma_Xgrav = sigma_x(N_common_ACC_Scale +1 : end ,1);
%-------------------------------------------------------------------------- 
else
Xgrav_neqreduced = 0;
sigma_Xgrav = 0;
end
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Normal Equations matrices assigned to the output grace_pod_struct structure array 
%-------------------------------------------------------------------------- 
grace_pod_struct.normal_equations_Nmatrix = NEQ_N; 
grace_pod_struct.normal_equations_umatrix = NEQ_u; 
grace_pod_struct.normal_equations_Nmatrix_reduced = NEQ_N_reduced; 
grace_pod_struct.normal_equations_umatrix_reduced = NEQ_u_reduced; 
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Special Cases of post-processing
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
end
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Normal Equations Combination Solution
%-------------------------------------------------------------------------- 
if test_grav_paramestim_yn == 0 || gravity_arc_sol_flag == 1
% cond_N_proscale_Fgrace_comb_neq = cond(NEQ_N)

% Combined NEQs solution
[Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xgrav_full, grace_pod_struct] = grace_comb_neqsol(grace_pod_struct, NEQ_N, NEQ_u, design_matrix_comb, reducedobservation_matrix, rangerate_empirical_struct, gravity_arc_sol_flag, Nglobal_scale); 

else
Xmatrix_orbit1 = 0; 
Xmatrix_orbit2 = 0; 
Xmatrix = 0; 
Xgrav_full = 0;     
end
%-------------------------------------------------------------------------- 

% Final Gravity solution
Xgrav = Xgrav_neqreduced;