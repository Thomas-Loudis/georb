function [orbit_model_struct] = orbit_model_init (main_config_fname,ic_satellites,src_version) 

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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Parallel Computing configuration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
param_keyword = 'parallel_computing_yn'; 
[parallel_computing_yn] = read_param_cfg(orbit_config_struct,param_keyword); 
 
param_keyword = 'num_cpu_cores'; 
[num_cpu_cores_char] = read_param_cfg(orbit_config_struct,param_keyword); 
num_cpu_cores = sscanf(num_cpu_cores_char,'%d %*'); 
 
orbit_model_struct.parallel_computing.parallel_computing_yn = parallel_computing_yn; 
orbit_model_struct.parallel_computing.num_cpu_cores         = num_cpu_cores; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
