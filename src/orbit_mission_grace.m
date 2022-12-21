function [out_dir_name] = orbit_mission_grace(main_config_fname, src_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_mission_grace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - config_filename     : Configuration file name
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            23 August  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB main settings :: Read main config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% % Satellite Data configuration file
% param_keyword = 'sat_data_filename';
% [sat_data_filename] = read_param_file(main_config_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode: GRACE missions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode: georb_mode :: 'orbit_mission' : Orbit determination of a single orbiting object (satellite, invidual orbiter)
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Satellite missions cases ::
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE or GRACE-FO mission' satellites
if test_mission_grace == 1
    satellite_1 = 'GRACE-A';
    satellite_2 = 'GRACE-B';
elseif test_mission_gracefo == 1
    satellite_1 = 'GRACE-C';
    satellite_2 = 'GRACE-D';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read IC configuration file to get the IC epochs list for the input orbiting mission/object
[ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_config_filename, satellite_1);
ic_data_satellite1 = ic_data_object1;

[ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_config_filename, satellite_2);
ic_data_satellite2 = ic_data_object1;

    
% Orbit analysis of GRACE series
[ic_n ic_m] =size(ic_data_satellite1);
for ic_i = 1 : ic_n
    estim_comb = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE-A / GRACE-C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current IC data
    ic_data_object_i = ic_data_satellite1(ic_i,:);
    
    % Orbit configuration structure
    [orbit_config_struct_GRACE1] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orbit Determination of GRACE-A/-C
    write_data = 1;
    [orbit_config_G1, sat1_orbit_matrix, sat1_orbit_rms, sat1_veqZ_matrix, sat1_veqP_matrix, sat1_OBS_matrix, sat1_Xparam_aposteriori] = orbit_object(orbit_config_struct_GRACE1, write_data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read config file for MJDo
    [orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_G1);    
    if test_mission_grace == 1 
        %OUT_foldername_GRACE = 'GRACE';
        OUT_foldername_GRACE = sprintf('%s%s%d','GRACE','_',fix(IC_MJDo));
    elseif test_mission_gracefo == 1
        %OUT_foldername_GRACE = 'GRACE-FO';    
        OUT_foldername_GRACE = sprintf('%s%s%d','GRACE-FO','_',fix(IC_MJDo));
    end    
    %[status, message, messageid] = rmdir(OUT_fname,'s');
    %[status, message, messageid] = mkdir(OUT_fname);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Move output files to GRACE/GRACE-FO folder name 
    if estim_comb == 0
    write_results_dir(orbit_config_G1, OUT_foldername_GRACE);
    end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE-B / GRACE-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRACE paired satellite :: Test of common MJD time 

    % Current IC data
    ic_data_object_i = ic_data_satellite2(ic_i,:);    

    % Orbit configuration structure
    [orbit_config_struct_GRACE2] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orbit Determination of GRACE-B/-D
    write_data = 1;
    [orbit_config_G2, sat2_orbit_matrix, sat2_orbit_rms, sat2_veqZ_matrix, sat2_veqP_matrix, sat2_OBS_matrix, sat2_Xparam_aposteriori] = orbit_object(orbit_config_struct_GRACE2, write_data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % Move output files to GRACE/GRACE-FO folder name 
    if estim_comb == 0
    write_results_dir(orbit_config_G2, OUT_foldername_GRACE);
    end
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
% KBR/LRI data residuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % K-Band Ranging (KBR) data analysis 
    intersat_obs = 'intersat_KBR';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, KBR_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs); 
    KBR_rms_res_rangerate = KBR_rms_resrangerate;
    KBR_rms_res_range     = rms_resrange;
    KBR_range_residuals     = resrange;
    KBR_rangerate_residuals = resrangerate;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LRI data analysis      
    intersat_obs = 'intersat_LRI';
    [biasrange, rangerate, rangeaccl, KBRbias, nonbiasrange, ...
          resrange, resrangerate, dresrange, dresrangerate, ...
          rms_resrange, LRI_rms_resrangerate, rms_dresrange, rms_dresrangerate]...
          = grace_kbr_analysis(main_config_fname, orbit_config_G1, orbit_config_G2, orbcGA, orbcGB, intersat_obs);      
    LRI_rms_res_rangerate = LRI_rms_resrangerate;
    LRI_rms_res_range     = rms_resrange;
    LRI_range_residuals     = resrange;
    LRI_rangerate_residuals = resrangerate;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1>0
fprintf('%s \n', 'Intersatellite Ranging residuals:');
fprintf('%s%21.6f', 'LRI range residuals: RMS (mm):          ',LRI_rms_res_range * 10^3);
fprintf('\n');
fprintf('%s%17.6f', 'LRI range-rate residuals: RMS (μm/sec): ',LRI_rms_res_rangerate * 10^6);
fprintf('\n');
fprintf('%s%17.6f', 'KBR range-rate residuals: RMS (μm/sec): ',KBR_rms_res_rangerate * 10^6);
fprintf('\n');
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Intersatellite Observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_matrix = LRI_rangerate_residuals;
data_functional = 'LRI range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);

data_matrix = LRI_range_residuals;
data_functional = 'LRI range residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);

data_matrix = KBR_rangerate_residuals;
data_functional = 'KBR range-rate residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct_GRACE1, data_matrix, data_functional, reference_frame);
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
    rms_lri(1,2) = LRI_rms_resrangerate;
    [fid] = write_georb_statistics(orbit_config_G1, orbit_config_G2, rms_obs, rms_extorb, rms_kbr, rms_lri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save GRACE orbits and KBR files to one archive
    [status,message,messageid] = movefile('*.out',OUT_foldername_GRACE);
    [status,message,messageid] = movefile('*.orb',OUT_foldername_GRACE);
    %[status,message,messageid] = movefile(out_dir_name_1,OUT_foldername_GRACE);
    %[status,message,messageid] = movefile(out_dir_name_2,OUT_foldername_GRACE);
    out_dir_name = OUT_foldername_GRACE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temporary overall output directory
    if ic_i == 1
        OUT_foldername_results_inprogress = 'results_in-progress';
        [status, message, messageid] = rmdir(OUT_foldername_results_inprogress);
        [status, message, messageid] = mkdir(OUT_foldername_results_inprogress);
    end
    [status,message,messageid] = movefile('*.out',OUT_foldername_results_inprogress,'f');
    [status,message,messageid] = movefile(out_dir_name,OUT_foldername_results_inprogress,'f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

out_dir_name = OUT_foldername_results_inprogress;    

