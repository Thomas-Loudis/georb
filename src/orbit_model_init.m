function [orbit_model_struct] = orbit_model_init (main_config_fname,ic_data_matrix,src_version)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model_init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit model intialisation; Create orbit/forces matrix for first epoch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               7 April 2025
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ic_data_orbiter1 = ic_data_matrix(1).ic_data;
ic_data_orbiter1_epoch0 = ic_data_orbiter1(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit configuration structure
[orbit_config_struct] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_orbiter1_epoch0, src_version);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit modelling : common models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to_pod = tic;
% Orbit Model structure matrix: Forces model, methods, data, parameters
orbit_model_struct.src_version = src_version;
[orbit_model_struct] = orbit_model (orbit_config_struct,orbit_model_struct);
orbit_model_struct.orbit_config_struct = orbit_config_struct;
orbit_model_struct.georb_config_fname = main_config_fname;
% fprintf('%s %.3f \n', 'Time (min):  call orbit_model :',toc(to_pod)/60);
fprintf('%s \n', 'Orbit Model initialisation :: Force model, Methods settings, Data preprocessing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
