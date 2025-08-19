function [orbit_model_struct] = orbit_model_init (main_config_fname,ic_satellites,src_version)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model_init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit model intialisation; Create orbit/forces matrix for first epoch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_model_struct  : Orbit model structure array 
% - ic_satellites       : Initial Conditions for all satellites for all days in structure array
%                       Structure Format: ic_satellites.epochs.ic_data
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               7 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read main configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit modelling configuration file
param_keyword = 'orb_config_filename';
[orbit_model_filename] = read_param_file(main_config_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ic_data_orbiter1 = ic_satellites(1).ic_data;
ic_data_orbiter1 = ic_satellites(1).epochs.ic_data;
ic_data_orbiter1_epoch0 = ic_data_orbiter1(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit configuration structure
[orbit_config_struct] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_orbiter1_epoch0, src_version);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main structure array for Orbit modelling : orbit_model_struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to_pod = tic;
% Orbit Model structure matrix: Forces model, methods, data, parameters
orbit_model_struct.src_version = src_version;
[orbit_model_struct] = orbit_model (orbit_config_struct,orbit_model_struct);
orbit_model_struct.orbit_config = orbit_config_struct;
orbit_model_struct.georb_config_fname = main_config_fname;
% fprintf('%s %.3f \n', 'Time (min):  call orbit_model :',toc(to_pod)/60);
fprintf('%s \n', 'Orbit Model initialisation :: Force model, Methods settings, Data preprocessing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
