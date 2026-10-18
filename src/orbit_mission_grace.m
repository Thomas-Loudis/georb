function [out_dir_name] = orbit_mission_grace(orbit_model_struct, ic_data_matrix) 

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
%  GRACE missions orbit determination and gravity field recovery  
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
% 23/05/2023, Thomas Loudis Papanikolaou 
%             Upgrade to support the gravity field parameters solution 
% 07/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
% 28/09/2026  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%--------------------------------------------------------------------------  
% GEORB main settings :: Read main config file 
%--------------------------------------------------------------------------  
main_config_fname = orbit_model_struct.georb_config_fname; 
 
% GEORB mode  
param_keyword = 'georb_mode'; 
[georb_mode] = read_param_file(main_config_fname,param_keyword); 
 
% Name ID of "satellite mission", "satellite" or "orbiting object"   
param_keyword = 'orbiting_object_name'; 
[orbiting_object_name] = read_param_file(main_config_fname,param_keyword); 
 
% Orbit modelling configuration files 
param_keyword = 'orb_config_filename'; 
[orbit_model_filename] = read_param_file(main_config_fname,param_keyword); 
 
src_version = orbit_model_struct.src_version; 
%--------------------------------------------------------------------------   
 
%--------------------------------------------------------------------------   
% Combined Parameter estimation settings  
%--------------------------------------------------------------------------   
observation_model_matrix = orbit_model_struct.observation_model_matrix; 

COMBESTIM_combparamestim_01 = observation_model_matrix.COMBESTIM_combparamestim_01; 
 
% Number of estimation occurences (If 0 :: Combined-Estimator not applied)    
COMBESTIM_Nestim_comb = observation_model_matrix.COMBESTIM_Nestim_comb; 
 
% Intersatellite Observations :: Laser (LRI) or K-band (KBR) 
COMBESTIM_intersat_obs = observation_model_matrix.COMBESTIM_intersat_obs;  
%--------------------------------------------------------------------------   
 
%--------------------------------------------------------------------------  
% Gravity Field parameter estimation approach  
%--------------------------------------------------------------------------  
% Gravity Field structure array 
gravity_field_struct = orbit_model_struct.gravity_field;
% Gravity solution y/n
gravity_solution_yn = gravity_field_struct.gravity_solution_yn; 
test_status_yn = strcmp(gravity_solution_yn,'y');
COMBESTIM_gravity_param = 0;
if test_status_yn == 1    
COMBESTIM_gravity_param = 1;
end     
%-------------------------------------------------------------------------- 
% Gravity Solution based on current orbit arc 
% to be applied for special cases only (simulation, daily solutions, ...)
%-------------------------------------------------------------------------- 
% 0. Arc-based gravity field solution is not computed and not written to ouput file; 
%    NEQs matrices are computed only and provided to the output
% 1. Arc-based gravity field solution is computed and written to ouput file
gravity_arc_sol_flag = gravity_field_struct.gravity_arc_sol_flag ;
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------  
% Accelerometers scale factors as global or arc parameters 
%--------------------------------------------------------------------------  
% Global parameters are included in the Normal Equations (NEQs) matrices
% following the pre-elimination of the orbit arc parameters. 
% Applied only in the case of gravity parameter estimation.
acc_scale_global_yn = orbit_model_struct.gravity_field.acc_scale_global_yn; 
test_status_yn = strcmp(acc_scale_global_yn,'y');    
if test_status_yn == 1    
acc_scale_common_01 = 1;
else
acc_scale_common_01 = 0;    
end     
%--------------------------------------------------------------------------  
 
%--------------------------------------------------------------------------  
% Mode: Satellite missions  
%--------------------------------------------------------------------------  
% Mode: georb_mode :: 'orbit_mission' : Orbit determination of a single orbiting object (satellite, invidual orbiter) 
test_orbit_mission = strcmp(georb_mode,'orbit_mission'); 
% Satellite missions cases :: 
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission'); 
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission'); 
test_mission_GRACEC  = strcmp(orbiting_object_name,'GRACE_C_mission'); 
test_mission_NGGM    = strcmp(orbiting_object_name,'NGGM_mission'); 
% GRACE or GRACE-FO mission' satellites 
if test_mission_grace == 1 || test_mission_GRACEC == 1 || test_mission_NGGM == 1 
    COMBESTIM_intersat_obs = 'KBR'; 
end 
intersat_obs_flag = COMBESTIM_intersat_obs; 
%--------------------------------------------------------------------------  
 
 
%--------------------------------------------------------------------------  
% POD of GRACE series 
%--------------------------------------------------------------------------   
% a-priori Initial Conditions (IC) of GRACE satellites 
ic_data_satellite1 = ic_data_matrix(1).ic_data; 
ic_data_satellite2 = ic_data_matrix(2).ic_data; 
[ic_n, ic_m] = size(ic_data_satellite1); 
ic_i = 1;

%-------------------------------------------------------------------------- 
% GRACE-A / GRACE-C orbit model array
%-------------------------------------------------------------------------- 
% Current IC data
ic_data_object_i = ic_data_satellite1(ic_i,:);
% Orbit configuration structure
[orbit_config_struct_GRACE1] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
% config_file = 'orbit_config_struct_GRACE1.txt'
% [config_file] = write_config_struct2file(orbit_config_struct_GRACE1,config_file);
% Orbit Model
[orbit_model_struct] = orbit_model (orbit_config_struct_GRACE1,orbit_model_struct);
% Orbit Data reading and preprocessing per satellite per date
[orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, orbit_config_struct_GRACE1, src_version, orbit_model_struct);
orbit_model_matrix_GRACE1 = orbit_model_struct;
orbit_model_matrix_GRACE1.orbit_config = orbit_config_struct_GRACE1;

%-------------------------------------------------------------------------- 
% GRACE-B / GRACE-D orbit model array
%-------------------------------------------------------------------------- 
% GRACE paired satellite
% Current IC data
ic_data_object_i = ic_data_satellite2(ic_i,:);
% Orbit configuration structure
[orbit_config_struct_GRACE2] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
% config_file = 'orbit_config_struct_GRACE2.txt'
% [config_file] = write_config_struct2file(orbit_config_struct_GRACE2,config_file);
% Orbit modelling : IC update
[orbit_model_struct] = orbit_model_ic (orbit_config_struct_GRACE2, orbit_model_struct);
% Orbit Data reading and preprocessing per satellite per date
[orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, orbit_config_struct_GRACE2, src_version, orbit_model_struct);
orbit_model_matrix_GRACE2 = orbit_model_struct;
orbit_model_matrix_GRACE2.orbit_config = orbit_config_struct_GRACE2;

%-------------------------------------------------------------------------- 
% grace_pod_struct initialisation
grace_pod_struct.grace1_pod.orbit_model = orbit_model_matrix_GRACE1;
grace_pod_struct.grace2_pod.orbit_model = orbit_model_matrix_GRACE2;

% orbit_matrices initialisation
[orbit_matrices] = orbit_matrices_init();

grace_pod_struct.grace1_pod.orbit_matrices = orbit_matrices;
grace_pod_struct.grace2_pod.orbit_matrices = orbit_matrices;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% POD step 1 :: individual orbit determination; non-combined parameter estimation
%-------------------------------------------------------------------------- 
write_data = 0;
intersat_ranging_flag = 'y';
[grace_pod_struct, out_dir_name] = grace_pod(grace_pod_struct, intersat_ranging_flag, intersat_obs_flag, write_data); 

% Write computed GRACE data to directory    
neq_flag = 0; 
[OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag, neq_flag); 
out_dir_name = OUT_fname_mission_mjd; 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Move apriori individual POD results to apriori directory  
OUT_fname_mission_mjd = out_dir_name; 
if COMBESTIM_combparamestim_01 == 1 
% folder_dest = 'Orbits_pod_apriori' 
OUT_dir_step1 = sprintf('%s%s',OUT_fname_mission_mjd,'_pod_step1'); 
[status,message,messageid] = movefile(OUT_fname_mission_mjd,OUT_dir_step1); 
% [status,message,messageid] = movefile('GRACE*',folder_dest); 
[status,message,messageid] = movefile('*.out',OUT_dir_step1); 
[status,message,messageid] = movefile('*.orb',OUT_dir_step1); 
[status,message,messageid] = movefile('*.prm',OUT_dir_step1); 
POD_apriori_orbits_folder = OUT_dir_step1; 
end 
%-------------------------------------------------------------------------- 


%--------------------------------------------------------------------------  
% Combined Parameter estimation based on combination of orbit and
% inter-satellite ranging observations (range-rate, range) 
%--------------------------------------------------------------------------  
if COMBESTIM_combparamestim_01 == 1 
fprintf('\n'); 
fprintf('%s \n', 'Combined Parameter Estimation algorithm :: in-progress'); 
fprintf('%s \n', 'Combination is based on orbits and inter-satellite ranging observations'); 
fprintf('\n'); 

%-------------------------------------------------------------------------- 
% Gravity parameter estimation 
%-------------------------------------------------------------------------- 
% Intermediate steps flag for gravity parameter estimation (y/n)
grav_param_01 = 0;

%-------------------------------------------------------------------------- 
% Gravity step model 1 or 2
COMBESTIM_gravity_step_model_ab = 1; 
%-------------------------------------------------------------------------- 

Nestim_comb = COMBESTIM_Nestim_comb; 
if COMBESTIM_gravity_step_model_ab == 2 
Nestim_comb = Nestim_comb + 1; 
end 
 
%-------------------------------------------------------------------------- 
% Range-Rate Empirical Parameters :: Initialisation of matrices
[rangerate_empirical_struct] = prm_intersat_empirical(orbit_model_struct);
% Structure array update :: Range-Rate Empirical Parameters  
grace_pod_struct.rangerate_empirical_struct = rangerate_empirical_struct; 
%-------------------------------------------------------------------------- 

% Iterations Loop
for i_iter_estim = 1 : Nestim_comb 
 
%-------------------------------------------------------------------------- 
% Update Gravity field and orbit model/configuration structure arrays      
%-------------------------------------------------------------------------- 
if COMBESTIM_gravity_param == 1  
if (COMBESTIM_gravity_step_model_ab == 1 && i_iter_estim >= 1) || (COMBESTIM_gravity_step_model_ab == 2 && i_iter_estim >= 2) 

fprintf('\n'); 
fprintf('%s \n', 'Gravity Partials :: computation in-progress'); 
fprintf('%s \n', 'Gravity partials are computed based on integration of variational equations'); 
fprintf('\n'); 
    
% Set orbit mode to orbit integration solution for orbit arc and partials    
orbit_pod_mode = 'orbit_propagation_veq'; 
 
% Initial Conditions (IC) values
% IC aposteriori update ::  
% Set to '0' to update IC 
% Set to '1' to add as IC corrections 
% Set to '-1' to avoid update IC via this call of grace_model_cor 
Xmatrix_flag = 0;      
IC_CRF_sat1 = grace_pod_struct.grace1_pod.orbit_model.IC_CRF; 
IC_CRF_sat2 = grace_pod_struct.grace2_pod.orbit_model.IC_CRF; 
sat1_Xparam_aposteriori = IC_CRF_sat1(1,2:end)'; 
sat2_Xparam_aposteriori = IC_CRF_sat2(1,2:end)'; 
  
% Modify configuration variables for Gravity parameters estimation
grav_param_01 = COMBESTIM_gravity_param; 
 
% Modify model :: Gravity Field structure array and orbit model/configuration arrays   
[grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_param_01, sat1_Xparam_aposteriori, sat2_Xparam_aposteriori, Xmatrix_flag); 
 
%-------------------------------------------------------------------------- 
% GRACE Partials / orbit propagation 
write_data = 0;
% InterSatellite Ranging Residuals computation based on initial observations (without empirical parameters) :: y/n
intersat_ranging_flag = 'n';
[grace_pod_struct, out_dir_name] = grace_pod(grace_pod_struct, intersat_ranging_flag, intersat_obs_flag, write_data);  
%-------------------------------------------------------------------------- 
end 
end 
% End of orbit Model/configuration Update and orbit partials computation  
%--------------------------------------------------------------------------  

%--------------------------------------------------------------------------  
% Combined Estimator algorithm :: Normal Equations Combination 
%--------------------------------------------------------------------------  
% NEQ_write 1/0 : Write/Exclude individual Normal Equations matrices per Observation set to ascii files
NEQ_write = 0;

% NEQs Combination and Solution
[NEQ_N, NEQ_u, NEQ_N_reduced, NEQ_u_reduced, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xgrav, sigma_Xgrav, grace_pod_struct] = grace_comb_neq(grace_pod_struct, acc_scale_common_01, NEQ_write, gravity_arc_sol_flag);
%--------------------------------------------------------------------------  

%-------------------------------------------------------------------------- 
if NEQ_write > 0 
OUT_foldername_ESTIM = sprintf('%s%d','neq_matrices_', i_iter_estim); 
[status,message,messageid] = mkdir(OUT_foldername_ESTIM); 
[status,message,messageid] = movefile('*.est',OUT_foldername_ESTIM); 
[status,message,messageid] = movefile('*.neq',OUT_foldername_ESTIM); 
end 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Update orbit model/configuration parameters for orbit integration iteration 
%-------------------------------------------------------------------------- 
if grav_param_01 == 0 || (grav_param_01 == 1 && gravity_arc_sol_flag == 1)
% Set mode of orbital equations integration solution  
if i_iter_estim == Nestim_comb 
% orbital integration solution of Equation of motion (orbtit only)    
orbit_pod_mode = 'orbit_propagation_eqm'; 
else 
% orbital integration solution of Equation of motion and Variational Equations (orbit and partials)    
orbit_pod_mode = 'orbit_propagation_veq'; 
end 
% IC Correction values update 
Xmatrix_flag = 1; 
% orbit model/configuration arrays update  
grav_update_flag = 0; 
[grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_update_flag, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix_flag);
end 
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Gravity Field parameters :: gravity model update
% Applied only for special cases (simulations, daily graity solutions)
%-------------------------------------------------------------------------- 
% Gravity Field solution update is handled individually via function gravity_param_aposteriori 
if grav_param_01 == 1 && gravity_arc_sol_flag == 1
fprintf('\n'); 
fprintf('%s \n', 'Gravity field solution :: Write the arc-based solution to output '); 
fprintf('\n'); 

% - write_gravity_flag  : Gravity sperical harmoncis coefficients to be written
% 1. Estimated Signal (Estimated coefficients corrections)
% 2. Final spherical harmonics coefficients :: Estimated Signal + Reference gravity model
write_gravity_flag = 2;
orbit_model_matrix = orbit_model_matrix_GRACE1;

% Gravity parameters estimated corrections to harmonic coefficients
[orbit_model_matrix, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori(Xgrav, sigma_Xgrav, orbit_model_matrix, write_gravity_flag); 

% Update Gravity Field model background for GRACE orbit model struct 
orbit_model_matrix_GRACE1.gravity_field = orbit_model_matrix.gravity_field; 
orbit_model_matrix_GRACE2.gravity_field = orbit_model_matrix.gravity_field; 

% grace_pod strucutre array update  
grace_pod_struct.grace1_pod.orbit_model = orbit_model_matrix_GRACE1; 
grace_pod_struct.grace2_pod.orbit_model = orbit_model_matrix_GRACE2; 
end
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% GRACE orbit propagation and partials 
if COMBESTIM_gravity_param == 0 || (grav_param_01 == 1 && gravity_arc_sol_flag == 1 && i_iter_estim == Nestim_comb)
write_data = 0;
intersat_ranging_flag = 'n';
[grace_pod_struct, out_dir_name] = grace_pod(grace_pod_struct, intersat_ranging_flag, intersat_obs_flag, write_data);   
end  
%-------------------------------------------------------------------------- 

end 
% End of Combined parameter estimator iterations 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Write computed GRACE data to directory in georb format 
neq_flag = 1; 
[OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag, neq_flag); 
if COMBESTIM_combparamestim_01 == 1 
[status,message,messageid] = movefile(POD_apriori_orbits_folder,OUT_fname_mission_mjd); 
[status,message,messageid] = movefile('neq_matrices_*',OUT_fname_mission_mjd);
if grav_param_01 == 1 && gravity_arc_sol_flag == 1
% [status,message,messageid] = movefile('*.gfc',OUT_foldername_ESTIM); 
[status,message,messageid] = movefile(gravity_model_filename,OUT_fname_mission_mjd); 
end 
end 
out_dir_name = OUT_fname_mission_mjd; 
%-------------------------------------------------------------------------- 
 
end 
% END of Combined Estimator algorithm 
%-------------------------------------------------------------------------- 
 

%--------------------------------------------------------------------------
% In-progress data output directory 
results_dir_name = 'results_in-progress'; 
results_dir_name_isfolder = isfolder(results_dir_name); 
if results_dir_name_isfolder == 0 
[status, message, messageid] = mkdir(results_dir_name); 
end 
% Move written output folders/files to results directory 
results_dir_path = fullfile(pwd,results_dir_name); 
[status,message,messageid] = movefile(out_dir_name,results_dir_path); 
%-------------------------------------------------------------------------- 
% end 
out_dir_name = results_dir_name;     
