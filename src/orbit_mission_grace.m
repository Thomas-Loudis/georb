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

% IC configuration file
param_keyword = 'ic_config_filename';
[ic_config_filename] = read_param_file(main_config_fname,param_keyword);

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

% Combined solution approach to Gravity and Orbit parameters estimation:  
COMBESTIM_comb_solution = observation_model_matrix.COMBESTIM_comb_solution;

% Combined solution Observations:  
COMBESTIM_OBS = observation_model_matrix.COMBESTIM_OBS;

% Intersatellite Observations :: Laser (LRI) or K-band (KBR)
COMBESTIM_intersat_obs_LRI = observation_model_matrix.COMBESTIM_intersat_obs_LRI;

% Weighted Estimation solution approaches
COMBESTIM_weight = observation_model_matrix.COMBESTIM_weight;

% Integration Step size reduced for LRI/KBR observations combined estimation
COMBESTIM_integration_step = observation_model_matrix.COMBESTIM_integration_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation approach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Status :: Not supported by current version 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMBESTIM_gravity_step = gravity_field_determination_matrix.COMBESTIM_gravity_step;
% Gravity parameters partials (in case COMBESTIM_gravity_step = 0)
COMBESTIM_gravity_partials = gravity_field_determination_matrix.COMBESTIM_gravity_partials;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode: GRACE missions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode: georb_mode :: 'orbit_mission' : Orbit determination of a single orbiting object (satellite, invidual orbiter)
test_orbit_mission = strcmp(georb_mode,'orbit_mission');
% Satellite missions cases ::
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission');
% GRACE or GRACE-FO mission' satellites
if test_mission_grace == 1
    satellite_1 = 'GRACE-A';
    satellite_2 = 'GRACE-B';
elseif test_mission_gracefo == 1
    satellite_1 = 'GRACE-C';
    satellite_2 = 'GRACE-D';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD of GRACE series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a-priori Initial Conditions of GRACE satellites
ic_data_satellite1 = ic_data_matrix(1).ic_data;
ic_data_satellite2 = ic_data_matrix(2).ic_data;
[ic_n ic_m] = size(ic_data_satellite1);
% IC loop start
for ic_i = 1 : ic_n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD (step 1 :: individually; non-combined estimation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_data = 1;
    if COMBESTIM_combparamestim_01 == 0
        write_data = 1;        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE-A / GRACE-C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current IC data
    ic_data_object_i = ic_data_satellite1(ic_i,:);  
    % Orbit configuration structure
    [orbit_config_struct_GRACE1] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version); 
    % Orbit Model 
    [orbit_model_struct] = orbit_model (orbit_config_struct_GRACE1,orbit_model_struct);
    % Orbit Data reading and preprocessing per satellite per date
    [orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, orbit_config_struct_GRACE1, src_version, orbit_model_struct);
    orbit_model_matrix_GRACE1 = orbit_model_struct;
    % Orbit Determination of GRACE-A/-C
    [orbit_config_G1, sat1_orbit_matrix, sat1_orbit_rms, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, sat1_Xparam_aposteriori, sat1_OBS_residuals, sat1_extorb] = orbit_object(orbit_config_struct_GRACE1, write_data, orbit_model_struct);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE-B / GRACE-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE paired satellite 
    % Current IC data
    ic_data_object_i = ic_data_satellite2(ic_i,:);    
    % Orbit configuration structure
    [orbit_config_struct_GRACE2] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
    % Orbit modelling : IC update
    [orbit_model_struct] = orbit_model_ic (orbit_config_struct_GRACE2, orbit_model_struct);
    % Orbit Data reading and preprocessing per satellite per date
    [orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, orbit_config_struct_GRACE2, src_version, orbit_model_struct);
    orbit_model_matrix_GRACE2 = orbit_model_struct;
    % Orbit Determination of GRACE-B/-D
    [orbit_config_G2, sat2_orbit_matrix, sat2_orbit_rms, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, sat2_Xparam_aposteriori, sat2_OBS_residuals, sat2_extorb] = orbit_object(orbit_config_struct_GRACE2, write_data, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nepochs, Nelements, N3] = size(sat1_orbit_matrix);
orbcGA = sat1_orbit_matrix(:,:,1);
orbcGB = sat2_orbit_matrix(:,:,1);

rms_obs_GA    = sat1_orbit_rms(1,1:3);
rms_orbitalGA = sat1_orbit_rms(2,1:3);
rms_obs_GB    = sat2_orbit_rms(1,1:3);
rms_orbitalGB = sat2_orbit_rms(2,1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intersatellite Ranging KBR/LRI data residuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRI_obs_analysis_01 = COMBESTIM_intersat_obs_LRI;
if test_mission_grace == 1
    LRI_obs_analysis_01 = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % K-Band Ranging (KBR) data residuals 
    intersat_obs = 'intersat_KBR';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, KBR_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_struct); 
    KBR_rms_res_rangerate = KBR_rms_resrangerate;
    KBR_rms_res_range     = rms_resrange;
    KBR_range_residuals     = resrange;
    KBR_rangerate_residuals = resrangerate;
    KBR_intersat_observation_data = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LRI_obs_analysis_01 == 1
    % Laser Ranging Interferometry LRI data residuals      
    intersat_obs = 'intersat_LRI';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, LRI_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_struct);      
    LRI_rms_res_rangerate = LRI_rms_resrangerate;
    LRI_rms_res_range     = rms_resrange;
    LRI_range_residuals     = resrange;
    LRI_rangerate_residuals = resrangerate;    
    LRI_intersat_observation_data = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)];  
else
    LRI_rms_res_rangerate = 0;
    LRI_rms_res_range     = 0;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Satellite Observations matrix   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LRI_obs_analysis_01 == 1
% LRI observations
intersat_observation_data = LRI_intersat_observation_data;
else
% KBR observations
intersat_observation_data = KBR_intersat_observation_data;
end
intersat_range_residuals = resrange;
intersat_rangerate_residuals = resrangerate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1>0
fprintf('%s \n', 'Intersatellite Ranging residuals:');
if LRI_obs_analysis_01 == 1
fprintf('%s%21.6f', 'LRI range residuals: RMS (mm):          ',LRI_rms_res_range * 10^3);
fprintf('\n');
fprintf('%s%17.6f', 'LRI range-rate residuals: RMS (μm/sec): ',LRI_rms_res_rangerate * 10^6);
fprintf('\n');
end
fprintf('%s%17.6f', 'KBR range residuals: RMS (mm):          ',KBR_rms_res_range * 10^3);
fprintf('\n');
fprintf('%s%17.6f', 'KBR range-rate residuals: RMS (μm/sec): ',KBR_rms_res_rangerate * 10^6);
fprintf('\n');
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Intersatellite Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mission Directory: GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis
[OUT_fname_object_mjd, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_matrix_GRACE1);
[status, message, messageid] = rmdir(OUT_fname_mission_mjd,'s');
[status, message, messageid] = mkdir(OUT_fname_mission_mjd);

if LRI_obs_analysis_01 == 1
data_matrix = LRI_rangerate_residuals;
data_functional = 'LRI range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = LRI_range_residuals;
data_functional = 'LRI range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
end

data_matrix = KBR_rangerate_residuals;
data_functional = 'KBR range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = KBR_range_residuals;
data_functional = 'KBR range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Statistics
    rms_extorb(1,1:3) = rms_orbitalGA;
    rms_extorb(2,1:3) = rms_orbitalGB;
    rms_obs(1,1:3)    = rms_obs_GA;
    rms_obs(2,1:3)    = rms_obs_GB;
    rms_kbr = zeros(1,3);
    rms_lri = zeros(1,3);
    rms_kbr(1,2) = KBR_rms_resrangerate;
if LRI_obs_analysis_01 == 1
    rms_lri(1,2) = LRI_rms_resrangerate;
end
[georb_data_name] = write_georb_statistics(orbit_config_G1, orbit_config_G2, rms_obs, rms_extorb, rms_kbr, rms_lri);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[OUT_fname_G1, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_matrix_GRACE1);
[OUT_fname_G2, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G2,orbit_model_matrix_GRACE2);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Move GRACE/GRACE-FO mission data files to one directory
    [status,message,messageid] = movefile('*.out',OUT_fname_mission_mjd);
    [status,message,messageid] = movefile('*.orb',OUT_fname_mission_mjd);
    [status,message,messageid] = movefile(OUT_fname_G1,OUT_fname_mission_mjd);
    [status,message,messageid] = movefile(OUT_fname_G2,OUT_fname_mission_mjd);
    out_dir_name = OUT_fname_mission_mjd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move apriori individual POD results to apriori directory 
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

% Update of configuration for observation combination    
if COMBESTIM_combparamestim_01 == 1 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration Update
% Set orbit mode to orbit integration solution for orbtit arc and partials   
orbit_pod_mode = 'orbit_propagation_veq';
% Gravity Field parameter estimation
grav_param_01 = COMBESTIM_gravity_param;
% GRACE-3 configuration update 
[orbit_config_struct_GRACE1,orbit_model_matrix_GRACE1] = config_struct_update(orbit_config_struct_GRACE1, orbit_pod_mode, sat1_Xparam_aposteriori, orbcGA, 0, grav_param_01, orbit_model_matrix_GRACE1);
% GRACE-4 configuration update
[orbit_config_struct_GRACE2,orbit_model_matrix_GRACE2] = config_struct_update(orbit_config_struct_GRACE2, orbit_pod_mode, sat2_Xparam_aposteriori, orbcGB, 0, grav_param_01, orbit_model_matrix_GRACE2);

if COMBESTIM_integration_step == 1
% Modify :: PARAM Stepsize :: VALUE 1 
param_keyword = 'Stepsize';
param_value = '1';
[orbit_config_struct_GRACE1] = write_configstruct_cor(orbit_config_struct_GRACE1, param_keyword, param_value);
[orbit_config_struct_GRACE2] = write_configstruct_cor(orbit_config_struct_GRACE2, param_keyword, param_value);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-A/C Partials
write_data = 0;
if COMBESTIM_Nestim_comb < 1
write_data = 1;
end
[orbit_config_G1, sat1_orbit_matrix, sat1_orbit_rms, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix_NOT, sat1_Xparam_aposteriori_NOT, sat1_OBS_residuals_NOT] = orbit_object(orbit_config_struct_GRACE1, write_data, orbit_model_matrix_GRACE1);

% GRACE-B/-D Partials
write_data = 0;
if COMBESTIM_Nestim_comb < 1
write_data = 1;
end
[orbit_config_G2, sat2_orbit_matrix, sat2_orbit_rms, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix_NOT, sat2_Xparam_aposteriori_NOT, sat2_OBS_residuals_NOT] = orbit_object(orbit_config_struct_GRACE2, write_data, orbit_model_matrix_GRACE2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated Orbit matrices of GRACE satellites 
orbcGA = sat1_orbit_matrix(:,:,1);
orbcGB = sat2_orbit_matrix(:,:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR/LRI data residuals update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % K-Band Ranging (KBR) data residuals 
    intersat_obs = 'intersat_KBR';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, KBR_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_matrix_GRACE1); 
    KBR_rms_res_rangerate = KBR_rms_resrangerate;
    KBR_rms_res_range     = rms_resrange;
    KBR_intersat_observation_data = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)];    

if LRI_obs_analysis_01 == 1
    % Laser Ranging Interferometry (LRI) data residuals      
    intersat_obs = 'intersat_LRI';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, LRI_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_matrix_GRACE1);      
    LRI_rms_res_rangerate = LRI_rms_resrangerate;
    LRI_rms_res_range     = rms_resrange;
    LRI_range_residuals     = resrange;
    LRI_rangerate_residuals = resrangerate;        
    LRI_intersat_observation_data = [rangerate(:,1) nonbiasrange(:,2) rangerate(:,2) rangeaccl(:,2)];  
else
    LRI_rms_res_rangerate = 0;
    LRI_rms_res_range     = 0;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Satellite Observations matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LRI_obs_analysis_01 == 1
% LRI observations
intersat_observation_data = LRI_intersat_observation_data;
else
% KBR observations
intersat_observation_data = KBR_intersat_observation_data;
end    
intersat_range_residuals = resrange;
intersat_rangerate_residuals = resrangerate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% End of Update of configuration and orbit partials computation 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations Weights 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted estimation solution approaches
weight_sol_opt = COMBESTIM_weight;
    
% Weights based on errors / residuals from step-1 orbit parameter estimation  
Cv_sat1_obs = [sat1_OBS_residuals(:,2); sat1_OBS_residuals(:,3); sat1_OBS_residuals(:,4)];
Cv_sat2_obs = [sat2_OBS_residuals(:,2); sat2_OBS_residuals(:,3); sat2_OBS_residuals(:,4)];  
% LRI range and range-rate obseravations
Cv_LRI_range = intersat_range_residuals(:,2);
Cv_LRI_rangerate = intersat_rangerate_residuals(:,2);
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];

% Fixed weights to sigma defined values
if weight_sol_opt == 3
sigma_obsorb = 5 * 10^-2;
sigma_range = 5 * 10^-3;
sigma_rangerate = 1 * 10^-6;
    
[d1_obs d2_obs] = size(Cv_sat1_obs);
[d1_obs2 d2_obs2] = size(Cv_sat2_obs);
[d1_lri d2_lri] = size(Cv_LRI_range);

Cv_sat1_obs = sigma_obsorb + zeros(d1_obs,1);
Cv_sat2_obs = sigma_obsorb + zeros(d1_obs2,1);
Cv_LRI_range = sigma_range + zeros(d1_lri,1);
Cv_LRI_rangerate = sigma_rangerate + zeros(d1_lri,1);
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];

elseif weight_sol_opt == 1
% Identity matrix
Cv_sat1_obs = 1;
Cv_sat2_obs = 1;
Cv_LRI_range = 1;
Cv_LRI_rangerate = 1;
Cv_LRI = [Cv_LRI_range Cv_LRI_rangerate];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation :: Combined Estimator
% Combined Parameter Estimation for GRACE satellites based on orbit and Laser intersatellite observations LRI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nestim_comb = COMBESTIM_Nestim_comb;
for i_iter_estim = 1 : Nestim_comb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-4 matrices for combined estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite observatios w.r.t. GRACE-4
[Xmatrix_LRI_sat2, Amatrix_rangerate_sat2, Wmatrix_rangerate_sat2, Amatrix_range_sat2, Wmatrix_range_sat2, NEQ_range_sat2, NEQ_rangerate_sat2] = estimator_orbit_intersat(orbcGA,orbcGB, sat1_veqZ_matrix,sat1_veqP_matrix, sat2_veqZ_matrix,sat2_veqP_matrix, intersat_observation_data, Cv_LRI);
%Xmatrix_LRI_1 = Xmatrix_LRI_sat2(1,1)

% Orbit pseudo-observations
[Xmatrix_obsorb_grace2,Xmatrix_alt_obsorb_grace2,Wmatrix_obsorb_grace2,Amatrix_obsorb_grace2, Cx, Cv, NEQn_grace2, NEQu_grace2] = estimator_orbit (orbcGB, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, Cv_sat2_obs,1);
%Xmatrix_obsorb_1_grace2 = Xmatrix_obsorb_grace2(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-3 matrices for combined estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite observatios w.r.t. GRACE-4
[Xmatrix_LRI_sat1, Amatrix_rangerate_sat1, Wmatrix_rangerate_sat1, Amatrix_range_sat1, Wmatrix_range_sat1, NEQ_range_sat1, NEQ_rangerate_sat1] = estimator_orbit_intersat(orbcGB,orbcGA, sat2_veqZ_matrix,sat2_veqP_matrix, sat1_veqZ_matrix,sat1_veqP_matrix, intersat_observation_data, Cv_LRI);
%Xmatrix_LRI_1_sat1 = Xmatrix_LRI_sat1(1,1)

% Orbit pseudo-observations
[Xmatrix_obsorb_grace1,Xmatrix_alt_obsorb_grace1,Wmatrix_obsorb_grace1,Amatrix_obsorb_grace1, Cx, Cv, NEQn_grace1, NEQu_grace1] = estimator_orbit (orbcGA, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, Cv_sat1_obs,1);
%Xmatrix_obsorb_1_grace1 = Xmatrix_obsorb_grace1(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfm_struct_glob = orbit_model_struct.gravity_field;
gfm_struct = gfm_struct_glob;
Nparam_common = 0;
if COMBESTIM_gravity_step > 0
grav_paramestim_yn = gfm_struct.param_estim_yn;
test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');   
if test_grav_paramestim_yn == 1   
N_param_GRAV = gfm_struct.parameters_number;
Nparam_common = N_param_GRAV;    
end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions of design matrix and Number of parameters, observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d1 d2]= size(Amatrix_obsorb_grace1);
[d3 d4]= size(Amatrix_obsorb_grace2);
Nparam_grace1 = d2;
Nparam_grace2 = d4;
Nobs_sat1 = d1;
Nobs_sat2 = d3;
Nparam_COMB = Nparam_grace1 + Nparam_grace2 - Nparam_common;

[d5 d6]= size(Amatrix_rangerate_sat1);
[d7 d8]= size(Amatrix_rangerate_sat2);
Nobs_rangerate = d7;
Nparam_rangerate = d8;

[ d9 d10]= size(Amatrix_range_sat1);
[d11 d12]= size(Amatrix_range_sat2);
Nobs_range = d9;
Nparam_range = d10;

N_obs_comb = Nobs_sat1 + Nobs_sat2 + Nobs_rangerate + Nobs_range;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amatrix = zeros(N_obs_comb , Nparam_COMB);
col_orbit1_1 = 1;
col_orbit1_2 = Nparam_grace1-Nparam_common;
col_orbit2_1 = col_orbit1_2 + 1;
col_orbit2_2 = col_orbit1_2 + Nparam_grace2-Nparam_common;
col_common_1 = col_orbit2_2 + 1;
col_common_2 = col_orbit2_2 + Nparam_common;

% Selected parameters to be estimated :: All or Golbal parameters (common)
if (Nparam_common > 0) && (COMBESTIM_comb_solution == 4)
param_col_start = col_common_1;
param_col_end = col_common_2;
else
param_col_start = 1;
param_col_end = Nparam_COMB;
end

% Orbit 1
A_orbit1 = zeros(Nobs_sat1 , param_col_end - (param_col_start-1));
A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1; 
% A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1(:,1: col_orbit1_2);

% Orbit 2
A_orbit2 = zeros(Nobs_sat2 , param_col_end - (param_col_start-1));
A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2; 
% A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2(:,1:Nparam_grace2-Nparam_common);

% Range-rate
col_common_0 = Nparam_rangerate-Nparam_common;
A_rangerate(1 : Nobs_rangerate , col_orbit1_1 : col_orbit1_2) = Amatrix_rangerate_sat1(:,1:col_common_0);
A_rangerate(1 : Nobs_rangerate , col_orbit2_1 : col_orbit2_2) = Amatrix_rangerate_sat2(:,1:col_common_0);

% Range
col_common_0 = Nparam_range-Nparam_common;
A_range(1 : Nobs_range, col_orbit1_1 : col_orbit1_2) = Amatrix_range_sat1(:,1:col_common_0);
A_range(1 : Nobs_range, col_orbit2_1 : col_orbit2_2) = Amatrix_range_sat2(:,1:col_common_0);

% Reduced-Observations matrix
b_orbit1 = Wmatrix_obsorb_grace1;
b_orbit2 = Wmatrix_obsorb_grace2;
b_rangerate = Wmatrix_rangerate_sat2;
b_range = Wmatrix_range_sat2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Estimation Solution based on Least Squares method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Approach-3
% Individual NEQ matrices per block/type of observations
% Orbit 1
[Xestim_orbit1, NEQ_N_orbit1, NEQ_u_orbit1] = estimator_neq_sol(A_orbit1, b_orbit1, Cv_sat1_obs);
% Orbit 2
[Xestim_orbit2, NEQ_N_orbit2, NEQ_u_orbit2] = estimator_neq_sol(A_orbit2, b_orbit2, Cv_sat2_obs);
% Range-rate
[Xestim_rangerate, NEQ_N_rangerate, NEQ_u_rangerate] = estimator_neq_sol(A_rangerate, b_rangerate, Cv_LRI_rangerate);
% Range
[Xestim_range, NEQ_N_range, NEQ_u_range] = estimator_neq_sol(A_range, b_range, Cv_LRI_range);

% Combined solution Observations:  
% Normal Equations :: Sum matrices
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

[NEQ_d1 NEQ_d2] = size(NEQ_N);
[NEQ_d1 NEQ_d2] = size(NEQ_u);

% Least Squares solution
tol2 = 30;
Xmatrix3 = lsqminnorm(NEQ_N, NEQ_u, tol2);
Xmatrix = Xmatrix3;
Xmatrix_Zo = Xmatrix(1:6,1);
[sz1, sz2] = size(Xmatrix);
Nparam_COMB_test = Nparam_COMB - sz1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEQ combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i_iter_estim == COMBESTIM_Nestim_comb
if ic_i == 1
[param_n, param_m] = size(NEQ_N);
[param_u, param_1] = size(NEQ_u);
NEQ_N_sum = zeros(param_n, param_n);
NEQ_u_sum = zeros(param_u, 1);
end
NEQ_N_sum = NEQ_N_sum + NEQ_N; 
NEQ_u_sum = NEQ_u_sum + NEQ_u; 

% Least Squares solution
tol2 = 30;
Xmatrix_NEQ_sum = lsqminnorm(NEQ_N_sum, NEQ_u_sum, tol2);
Xmatrix_6_NEQ_sum = Xmatrix_NEQ_sum(1:6,1);
if ic_i == ic_n   
Xmatrix = Xmatrix_NEQ_sum; 
Xmatrix_6_NEQ_sum = Xmatrix(1:6,1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit and Gravity Field parameters estimated corrections
if (Nparam_common > 0) && (COMBESTIM_comb_solution == 4)
% Gravity parameters corrections    
Xmatrix_gravparam = Xmatrix;
% Set orbit parameters corrections to zero
Xmatrix_orbit1 = zeros(col_orbit1_2 , 1); 
% Xmatrix_orbit2 = zeros(col_orbit2_1 : col_orbit2_2 , 1);
Xmatrix_orbit2 = zeros(col_orbit1_2 , 1);
else
% Orbit Parameters estimated corrections    
Xmatrix_orbit1 = Xmatrix(1 : col_orbit1_2 , 1);
Xmatrix_orbit2 = Xmatrix(col_orbit2_1 : col_orbit2_2 , 1);   
% Gravity parameters estimated corrections    
Xmatrix_gravparam = Xmatrix(col_common_1 : col_common_2 , 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not supported by current Version 
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT_foldername_ESTIM = sprintf('%s%d','param_estim_', i_iter_estim);
% [status, message, messageid] = mkdir(OUT_foldername_ESTIM);
% [status,message,messageid] = movefile('*.est',OUT_foldername_ESTIM);
% [status,message,messageid] = movefile('*.gfc',OUT_foldername_ESTIM);
% [status,message,messageid] = movefile('*.neq',OUT_foldername_ESTIM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update configuration parameters for orbit integration iteration
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
% Gravity Field parameter update
if COMBESTIM_gravity_step > 0
grav_param_01 = COMBESTIM_gravity_param;
if i_iter_estim == Nestim_comb
if COMBESTIM_gravity_param == 0 && COMBESTIM_gravity_partials == 1     
grav_param_01 = COMBESTIM_gravity_partials;
orbit_pod_mode = 'orbit_propagation_veq'; 
% Update configuration file :: parameters to be estimated
paramestim_update_config_01 = 1;
end
end
if COMBESTIM_gravity_step == 2
COMBESTIM_gravity_param = 1 ;   
orbit_pod_mode = 'orbit_propagation_veq';   
grav_param_01 = 1;
% Update configuration file :: parameters to be estimated
paramestim_update_config_01 = 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i_iter_estim == Nestim_comb
% orbital integration solution of Equation of motion (orbtit only)   
orbit_pod_mode = 'orbit_propagation_eqm';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-3 config update and orbit integration (EQM + VEQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmatrix = Xmatrix_orbit1;
% Add Xaposteriori to functions: param_aposteriori, orbit_pod, orbit_object 
if i_iter_estim == 1
% Xapriori = Xaposteriori;
Xapriori = sat1_Xparam_aposteriori (1 : Nparam_grace1-Nparam_common, 1);
else
Xapriori = Xparam_apriori_sat1;
end
% Xaposteriori = Xapriori + Xmatrix (obs and lri)
Xaposteriori = Xapriori + Xmatrix;
Xparam_apriori_sat1 = Xaposteriori; 
% GRACE-3 configuration update
[orbit_config_struct_GRACE1,orbit_model_matrix_GRACE1] = config_struct_update(orbit_config_struct_GRACE1, orbit_pod_mode, Xaposteriori, orbcGA, 0, grav_param_01,orbit_model_matrix_GRACE1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-4 config update and orbit integration (EQM + VEQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmatrix = Xmatrix_orbit2;
% Add Xaposteriori to functions: param_aposteriori, orbit_pod, orbit_object 
if i_iter_estim == 1
% Xapriori = Xaposteriori;
Xapriori = sat2_Xparam_aposteriori (1 : Nparam_grace2-Nparam_common, 1);
else
Xapriori = Xparam_apriori_sat2;
end
% Xaposteriori = Xapriori + Xmatrix (obs and lri)
Xaposteriori = Xapriori + Xmatrix;
Xparam_apriori_sat2 = Xaposteriori;
% GRACE-4 configuration structure update
[orbit_config_struct_GRACE2,orbit_model_matrix_GRACE2] = config_struct_update(orbit_config_struct_GRACE2, orbit_pod_mode, Xaposteriori, orbcGB, 0, grav_param_01,orbit_model_matrix_GRACE2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration of GRACE-A/C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_data = 0;
if i_iter_estim == Nestim_comb
write_data = 1;
end
[orbit_config_G1, sat1_orbit_matrix, sat1_orbit_rms, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix_NOT, sat1_Xparam_aposteriori_NOT, sat1_OBS_residuals_NOT,ext_orb,forces_accel,orbit_model_matrix_GRACE1] = orbit_object(orbit_config_struct_GRACE1, write_data, orbit_model_matrix_GRACE1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration of GRACE-B/-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_data = 0;
if i_iter_estim == Nestim_comb
write_data = 1;
end
[orbit_config_G2, sat2_orbit_matrix, sat2_orbit_rms, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix_NOT, sat2_Xparam_aposteriori_NOT, sat2_OBS_residuals_NOT,ext_orb,forces_accel,orbit_model_matrix_GRACE2] = orbit_object(orbit_config_struct_GRACE2, write_data, orbit_model_matrix_GRACE2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated estimated Orbit matrices of GRACE satellites 
orbcGA = sat1_orbit_matrix(:,:,1);
orbcGB = sat2_orbit_matrix(:,:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KBR/LRI data residuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % K-Band Ranging (KBR) data analysis 
    intersat_obs = 'intersat_KBR';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, KBR_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_matrix_GRACE1); 
    KBR_rms_res_rangerate = KBR_rms_resrangerate;
    KBR_rms_res_range     = rms_resrange;
    KBR_range_residuals     = resrange;
    KBR_rangerate_residuals = resrangerate;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LRI_obs_analysis_01 == 1
    % LRI data analysis      
    intersat_obs = 'intersat_LRI';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, LRI_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs, orbit_model_matrix_GRACE1);      
    LRI_rms_res_rangerate = LRI_rms_resrangerate;
    LRI_rms_res_range     = rms_resrange;
    LRI_range_residuals     = resrange;
    LRI_rangerate_residuals = resrangerate;  
else
    LRI_rms_res_rangerate = 0;
    LRI_rms_res_range     = 0;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1>0
fprintf('%s \n', 'Intersatellite Ranging residuals:');
if LRI_obs_analysis_01 == 1
fprintf('%s%21.6f', 'LRI range residuals: RMS (mm):          ',LRI_rms_res_range * 10^3);
fprintf('\n');
fprintf('%s%17.6f', 'LRI range-rate residuals: RMS (μm/sec): ',LRI_rms_res_rangerate * 10^6);
fprintf('\n');
end
fprintf('%s%17.6f', 'KBR range residuals: RMS (mm):          ',KBR_rms_res_range * 10^3);
fprintf('\n');
fprintf('%s%17.6f', 'KBR range-rate residuals: RMS (μm/sec): ',KBR_rms_res_rangerate * 10^6);
fprintf('\n');
fprintf('\n')
end

end
% End of Combined parameter estimator iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Intersatellite Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mission Directory: GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis
[OUT_fname_G1, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G1,orbit_model_matrix_GRACE1);
[OUT_fname_G2, OUT_fname_mission_mjd] = write_results_dir(orbit_config_G2,orbit_model_matrix_GRACE2);
[status, message, messageid] = rmdir(OUT_fname_mission_mjd,'s');
[status, message, messageid] = mkdir(OUT_fname_mission_mjd);

if LRI_obs_analysis_01 == 1
data_matrix = LRI_rangerate_residuals;
data_functional = 'LRI range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = LRI_range_residuals;
data_functional = 'LRI range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
end

data_matrix = KBR_rangerate_residuals;
data_functional = 'KBR range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);

data_matrix = KBR_range_residuals;
data_functional = 'KBR range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Residuals results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rms_obs_GA    = sat1_orbit_rms(1,1:3);
rms_orbitalGA = sat1_orbit_rms(2,1:3);
% rms_obs_GB    = sat2_orbit_rms(1,1:3);
rms_orbitalGB = sat2_orbit_rms(2,1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Orbit Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE1
orbc = orbcGA;
orbce = sat1_OBS_matrix;
orbc_pos = orbc(:,1:4);
% Celestial Reference Frame (GCRS)
[dorbc,rms_orbc,orbc_common] = compstat(orbce,orbc_pos);
rms_obs_GA = rms_orbc;
OBS_residuals = dorbc;

data_matrix = OBS_residuals;
data_functional = 'Orbit Observations residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_G1);

% GRACE2
orbc = orbcGB;
orbce = sat2_OBS_matrix;
orbc_pos = orbc(:,1:4);
% Celestial Reference Frame (GCRS)
[dorbc,rms_orbc,orbc_common] = compstat(orbce,orbc_pos);
rms_obs_GB = rms_orbc;
OBS_residuals = dorbc;

data_matrix = OBS_residuals;
data_functional = 'Orbit Observations residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE2, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_G2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters (Vector matrix) :: Orbital Parameters
% GRACE-1/3
data_matrix = Xparam_apriori_sat1';
data_functional = 'Parameters';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_G1);
% GRACE-2/4
data_matrix = Xparam_apriori_sat2';
data_functional = 'Parameters';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE2, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_G2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Statistics
    rms_extorb(1,1:3) = rms_orbitalGA;
    rms_extorb(2,1:3) = rms_orbitalGB;
    rms_obs(1,1:3)    = rms_obs_GA;
    rms_obs(2,1:3)    = rms_obs_GB;
    rms_kbr = zeros(1,3);
    rms_lri = zeros(1,3);
    rms_kbr(1,2) = KBR_rms_resrangerate;
if LRI_obs_analysis_01 == 1
    rms_lri(1,2) = LRI_rms_resrangerate;
end
    [georb_data_name] = write_georb_statistics(orbit_config_G1, orbit_config_G2, rms_obs, rms_extorb, rms_kbr, rms_lri);
    [status,message,messageid] = movefile(georb_data_name, OUT_fname_mission_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if COMBESTIM_combparamestim_01 == 1
[status,message,messageid] = movefile(POD_apriori_orbits_folder,OUT_fname_mission_mjd);
[status,message,messageid] = movefile('param_estim_*',OUT_fname_mission_mjd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Move GRACE/GRACE-FO mission data files to one directory
    [status,message,messageid] = movefile('*.out',OUT_fname_mission_mjd);
    [status,message,messageid] = movefile('*.orb',OUT_fname_mission_mjd);
    [status,message,messageid] = movefile(OUT_fname_G1,OUT_fname_mission_mjd);
    [status,message,messageid] = movefile(OUT_fname_G2,OUT_fname_mission_mjd);
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
