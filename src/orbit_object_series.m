function [out_dir_name] = orbit_object_series(orbit_model_struct, ic_data_matrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_object_series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  orbit_object_series is used for the orbit determination of objects series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - main_config_fname     : GEORB master configuration file name
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             23 August 2022
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

ic_data = ic_data_matrix.ic_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to output files :: y/n
write_data = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit analysis of Objects series loop
[ic_n ic_m] = size(ic_data);
for ic_i = 1 : ic_n
    % Current IC data
    ic_data_object_i = ic_data(ic_i,:);
    
    % Orbit configuration structure
    [config_struct] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
        
    % Orbit Model matrix
    [orbit_model_struct] = orbit_model (config_struct,orbit_model_struct);
    % [orbit_model_struct] = orbit_model_ic (config_struct, orbit_model_struct);

    % Orbit Data reading and preprocessing per satellite per date
    [orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, config_struct, src_version, orbit_model_struct);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orbit Determination of satellite/object i
    [orbit_config_fname, orbit_matrix, orbit_rms, veqZ_matrix, veqP_matrix, OBS_matrix, Xparam_aposteriori, OBS_residuals, extorb] = orbit_object(config_struct, write_data, orbit_model_struct);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Move output data files to the object output data directory
    [OUT_fname_object_mjd, OUT_fname_mission_mjd] = write_results_dir(orbit_config_fname,orbit_model_struct);
    [status,message,messageid] = movefile('*.out',OUT_fname_object_mjd);
    [status,message,messageid] = movefile('*.orb',OUT_fname_object_mjd);

    % in-progress data output directory
    results_dir_name = 'results_in-progress';
    results_dir_name_isfolder = isfolder(results_dir_name);
    if results_dir_name_isfolder == 0
    [status, message, messageid] = mkdir(results_dir_name);
    end    
    % Move written output folders/files to results directory
    results_dir_path = fullfile(pwd,results_dir_name);
    [status,message,messageid] = movefile(OUT_fname_object_mjd,results_dir_path);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir_name = results_dir_name;    
