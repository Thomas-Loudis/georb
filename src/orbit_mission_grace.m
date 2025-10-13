function [out_dir_name] = orbit_mission_grace(orbit_model_struct, ic_data_matrix)

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
% 23/05/2023, Thomas Loudis Papanikolaou
%             Upgrade to support the gravity field parameters solution
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB main settings :: Read main config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Parameter estimation settings (TEMP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observation_model_matrix = orbit_model_struct.observation_model_matrix;
gravity_field_determination_matrix = orbit_model_struct.gravity_field_determination_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMBESTIM_combparamestim_01 = observation_model_matrix.COMBESTIM_combparamestim_01;

% Number of estimation occurences (If 0 :: Combined-Estimator not applied)   
COMBESTIM_Nestim_comb = observation_model_matrix.COMBESTIM_Nestim_comb;

% Combined solution Observations:  
COMBESTIM_OBS = observation_model_matrix.COMBESTIM_OBS;

% Intersatellite Observations :: Laser (LRI) or K-band (KBR)
COMBESTIM_intersat_obs = observation_model_matrix.COMBESTIM_intersat_obs;

% Weighted Estimation solution approaches
COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation approach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Status :: Not supported by current version 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMBESTIM_gravity_step = gravity_field_determination_matrix.COMBESTIM_gravity_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode: GRACE missions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD of GRACE series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a-priori Initial Conditions of GRACE satellites
ic_data_satellite1 = ic_data_matrix(1).ic_data;
ic_data_satellite2 = ic_data_matrix(2).ic_data;
[ic_n, ic_m] = size(ic_data_satellite1);
% IC loop start
for ic_i = 1 : ic_n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-A / GRACE-C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-B / GRACE-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD (step 1 :: individually; non-combined estimation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_data = 0;
[grace_pod_struct, grace1_pod, grace2_pod, intersat_pod, out_dir_name] = grace_pod(orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2, intersat_obs_flag, write_data);
grace_pod_struct_step1 = grace_pod_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write computed GRACE data to directory   
neq_flag = 0;
[OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag, neq_flag);
out_dir_name = OUT_fname_mission_mjd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Parameter estimation based on orbit and inter-satellite ranging observations (range-rate, range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if COMBESTIM_combparamestim_01 == 1

% Gravity parameter estimation step-wise approach :: Steps
if COMBESTIM_gravity_step == 2
    COMBESTIM_gravity_param = 0;
end
if COMBESTIM_gravity_step == 0
    COMBESTIM_gravity_param = 0;
end
if COMBESTIM_gravity_step == 1
    COMBESTIM_gravity_param = 1;
end

grav_param_01 = 0;

% Observations Residuals from POD step 1 
grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals = grace_pod_struct_step1.grace1_pod.orbit_matrices.observation_residuals; 
grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals = grace_pod_struct_step1.grace2_pod.orbit_matrices.observation_residuals; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation :: Combined Estimator
% Combined Parameter Estimation for GRACE satellites based on orbit observations (GPS-based kinematic positions) and intersatellite range-rate observations (KBR, LRI)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nestim_comb = COMBESTIM_Nestim_comb;

COMBESTIM_gravity_step_model_ab = 1;
if COMBESTIM_gravity_step_model_ab == 2
Nestim_comb = Nestim_comb + 1;
end

for i_iter_estim = 1 : Nestim_comb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update orbit model/configuration for gravity field parameters and observation combination solution    
if COMBESTIM_gravity_param == 1 % COMBESTIM_combparamestim_01 == 1
if COMBESTIM_gravity_step_model_ab == 1 || (COMBESTIM_gravity_step_model_ab == 2 && i_iter_estim == 2)

% Set orbit mode to orbit integration solution for orbit arc and partials   
orbit_pod_mode = 'orbit_propagation_veq';
% Gravity Field parameter estimation
grav_param_01 = 0;
% Initial Conditions (IC) values
% Xparam_aposteriori = grace_pod_struct.grace1_pod.orbit_model.IC_CRF
sat1_Xparam_aposteriori = grace_pod_struct.grace1_pod.orbit_matrices.param_aposteriori;
sat2_Xparam_aposteriori = grace_pod_struct.grace2_pod.orbit_matrices.param_aposteriori;
% IC aposteriori 
Xmatrix_flag = 0;
% orbit model/configuration arrays update 
% [grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_param_01, sat1_Xparam_aposteriori, sat2_Xparam_aposteriori, Xmatrix_flag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Gravity Field solution
grav_param_01 = COMBESTIM_gravity_param;
% IC aposteriori 
Xmatrix_flag = 0;   % Xmatrix_flag = -1;
% orbit model/configuration arrays update 
[grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_param_01, sat1_Xparam_aposteriori, sat2_Xparam_aposteriori, Xmatrix_flag);
% Update Gravity Field structure array
if grav_param_01 == 1 
grav_field_paramestim_yn = 'y';
[orbit_model_matrix_GRACE1] = orbit_model_gravity (grav_field_paramestim_yn,orbit_model_matrix_GRACE1);
[orbit_model_matrix_GRACE2] = orbit_model_gravity (grav_field_paramestim_yn,orbit_model_matrix_GRACE2);
grace_pod_struct.grace1_pod.orbit_model = orbit_model_matrix_GRACE1;
grace_pod_struct.grace2_pod.orbit_model = orbit_model_matrix_GRACE2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE Partials / orbit propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_data = 0;
[grace_pod_struct, grace1_pod, grace2_pod, intersat_pod, out_dir_name] = grace_pod(orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2, intersat_obs_flag, write_data); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
% End of orbit Model/configuration Update and orbit partials computation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (pseudo-)Observations to orbits (GPS-based kinematic positions) from POD step 1 matrices   
grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace1_pod.orbit_matrices.observation_matrix; 
grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace2_pod.orbit_matrices.observation_matrix; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Estimator algorithm :: Design matrix, Normal Equations     
design_matrix_scale_orbit = 0; 
design_matrix_param_orbit = 0; 
NEQ_write = 0;
[NEQ_N, NEQ_u, NEQ_N_reduced, NEQ_u_reduced, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix, Xcommon, Xcommon_NEQreduced, grace_pod_struct] = grace_comb_neq(grace_pod_struct, design_matrix_scale_orbit, design_matrix_param_orbit, NEQ_write);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update orbit model/configuration parameters for orbit integration iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set mode of orbital equations integration solution 
if i_iter_estim == Nestim_comb
% orbital integration solution of Equation of motion (orbtit only)   
orbit_pod_mode = 'orbit_propagation_eqm';
else
% orbital integration solution of Equation of motion and Variational Equations (orbtit and partials)   
orbit_pod_mode = 'orbit_propagation_veq';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC Correction values update
% if grav_param_01 == 0
Xmatrix_flag = 1;
% orbit model/configuration arrays update 
% [grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_param_01, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix_flag); 
grav_update_flag = 0;
[grace_pod_struct, orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2] = grace_model_cor(grace_pod_struct, orbit_pod_mode, grav_update_flag, Xmatrix_orbit1, Xmatrix_orbit2, Xmatrix_flag); 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters :: gravity model update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field solution update is handled individually via function gravity_param_aposteriori
if grav_param_01 == 1
ic_apriori_01 = 1;
% [orbit_model_matrix_GRACE1, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori(Xcommon, ic_apriori_01, orbit_model_matrix_GRACE1);
[orbit_model_matrix_GRACE1, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori(Xcommon_NEQreduced, ic_apriori_01, orbit_model_matrix_GRACE1);
% Update Gravity Field model background for GRACE 2/4
orbit_model_matrix_GRACE2.gravity_field = orbit_model_matrix_GRACE1.gravity_field;

% save Cmatrix_aposteriori.neq  Cmatrix_aposteriori -ASCII -double
% save Smatrix_aposteriori.neq  Smatrix_aposteriori -ASCII -double
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NEQ_write > 0
OUT_foldername_ESTIM = sprintf('%s%d','neq_matrices_', i_iter_estim);
[status, message, messageid] = mkdir(OUT_foldername_ESTIM);
[status,message,messageid] = movefile('*.est',OUT_foldername_ESTIM);
[status,message,messageid] = movefile('*.neq',OUT_foldername_ESTIM);
if grav_param_01 == 1
% [status,message,messageid] = movefile('*.gfc',OUT_foldername_ESTIM);
[status,message,messageid] = movefile(gravity_model_filename,OUT_foldername_ESTIM);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if grav_param_01 == 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE orbit propagation and partials
write_data = 0;
[grace_pod_struct, grace1_pod, grace2_pod, intersat_pod, out_dir_name] = grace_pod(orbit_model_matrix_GRACE1, orbit_model_matrix_GRACE2, intersat_obs_flag, write_data);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Parameters
grace_pod_struct.grace1_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit1;
grace_pod_struct.grace2_pod.orbit_matrices.param_aposteriori = Xmatrix_orbit2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation Residuals update 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (pseudo-)Observations to orbits (GPS-based kinematic positions) from POD step 1 matrices   
grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace1_pod.orbit_matrices.observation_matrix; 
grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace2_pod.orbit_matrices.observation_matrix; 

% Observation matrices
G1_observation_matrix = grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix; 
G2_observation_matrix = grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix; 

orbc_G1 = grace_pod_struct.grace1_pod.orbit_matrices.orbit_crf;
orbc_G2 = grace_pod_struct.grace2_pod.orbit_matrices.orbit_crf;

orbc_G1_xyz = orbc_G1(:,1:4);
orbc_G2_xyz = orbc_G2(:,1:4);

[dorbc,rms_orbc,orbc_common] = compstat(G1_observation_matrix,orbc_G1_xyz); 
obs_residuals_G1 = dorbc;
[dorbc,rms_orbc,orbc_common] = compstat(G2_observation_matrix,orbc_G2_xyz);  
obs_residuals_G2 = dorbc;

grace_pod_struct.grace1_pod.orbit_matrices.observation_residuals = obs_residuals_G1; 
grace_pod_struct.grace2_pod.orbit_matrices.observation_residuals = obs_residuals_G2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end
% End of Combined parameter estimator iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write computed GRACE Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (pseudo-)Observations to orbits (GPS-based kinematic positions) from POD step 1 matrices   
grace_pod_struct.grace1_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace1_pod.orbit_matrices.observation_matrix; 
grace_pod_struct.grace2_pod.orbit_matrices.observation_matrix = grace_pod_struct_step1.grace2_pod.orbit_matrices.observation_matrix; 

% Normal Equations matrices
grace_pod_struct.normal_equations_Nmatrix = NEQ_N;
grace_pod_struct.normal_equations_umatrix = NEQ_u;
grace_pod_struct.normal_equations_Nmatrix_reduced = NEQ_N_reduced;
grace_pod_struct.normal_equations_umatrix_reduced = NEQ_u_reduced;

% Write computed GRACE data to directory in georb format
neq_flag = 1;
[OUT_fname_mission_mjd,OUT_data_foldername_G1,OUT_data_foldername_G2, grace_pod_struct] = grace_writedata(grace_pod_struct, intersat_obs_flag, neq_flag);
if COMBESTIM_combparamestim_01 == 1
[status,message,messageid] = movefile(POD_apriori_orbits_folder,OUT_fname_mission_mjd);
[status,message,messageid] = movefile('neq_matrices_*',OUT_fname_mission_mjd);
end
out_dir_name = OUT_fname_mission_mjd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% END of Combined Estimator algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-progress data output directory
results_dir_name = 'results_in-progress';
results_dir_name_isfolder = isfolder(results_dir_name);
if results_dir_name_isfolder == 0
[status, message, messageid] = mkdir(results_dir_name);
end
% Move written output folders/files to results directory
results_dir_path = fullfile(pwd,results_dir_name);
[status,message,messageid] = movefile(out_dir_name,results_dir_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
out_dir_name = results_dir_name;    
