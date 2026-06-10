function [orbit_config_struct, orbit_model_struct, orbit_matrices, OUT_data_foldername] = orbit_object(orbit_config_struct, write_data, orbit_model_struct) 

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
% Function: orbit_object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  orbit_object performs the individual orbit analysis (precise orbit 
%  determination, orbit propagation, simulation) of a single orbiting 
%  object (satellite, individual object)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - main_config_fname     : GEORB master configuration file name 
% - ic_data_object        :  
% 
% Output arguments: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                             23 August 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 07/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit Determination mode 
param_keyword = 'orbit_mode'; 
[orbit_mode] = read_param_cfg(orbit_config_struct,param_keyword); 
test = strcmp(orbit_mode,'orbit_propagation_veq'); 
if test == 1 
    fprintf('\n%s\n','Mode: Orbit Propagation EQM and VEQ');     
end 
test = strcmp(orbit_mode,'orbit_propagation_eqm'); 
if test == 1 
    fprintf('\n%s\n','Mode: Orbit Propagation EQM');     
end 
test = strcmp(orbit_mode,'orbit_determination'); 
if test == 1 
    fprintf('\n%s\n','Mode: Orbit Determination');     
end 
 
% Satellite or Orbiting OBject name 
param_keyword = 'orbiting_object_name'; 
[orbiting_object_name] = read_param_cfg(orbit_config_struct,param_keyword); 
fprintf('%s%s\n','Orbiting Object: ', orbiting_object_name); 
 
% Initial Epoch 
param_keyword = 'Initial_Epoch'; 
[param_value, param_line] = read_param_cfg(orbit_config_struct,param_keyword); 
Initial_Epoch = param_line; 
fprintf('%s%s\n','Initial Epoch: ',Initial_Epoch); 
fprintf('\n'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit Determination :: Call POD function orbit_pod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
test = strcmp(orbit_mode,'orbit_determination'); 
if test == 1 
fprintf('%s\n', 'Orbit Determination :: in-progress : Iterations of Orbit Integration & Parameter Estimation');    
else 
fprintf('%s\n', 'Orbit Propagation :: in-progress');    
end    
% fprintf('%s\n','Models and Data preprocessing :: In-progress');  
 
% Call function orbit_pod 
[orbit_model_struct, orbit_matrices] = orbit_pod(orbit_config_struct, orbit_model_struct); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
rms_orbital = orbit_matrices.rms_orbital; 
rms_orbc    = orbit_matrices.rms_orbc; 
rms_obs     = orbit_matrices.rms_observation_residuals; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
test = strcmp(orbit_mode,'orbit_determination'); 
if test == 1 
    % fprintf('\n'); 
    fprintf('%s %11.6f %11.6f %11.6f', 'Orbit residuals: RMS(XYZ): ',rms_obs(1:3)); 
    fprintf('\n\n'); 
end 
 
param_keyword = 'external_orbit_comp'; 
[external_orbit_comp_yn] = read_param_cfg(orbit_config_struct,param_keyword); 
test = strcmp(external_orbit_comp_yn,'y'); 
if test == 1 
    fprintf('%s \n', 'External Orbit Comparison:'); 
    fprintf('%s ', 'Orbit residuals: RMS(RTN): ' ); 
    fprintf('%11.6f ', rms_orbital); 
    fprintf('\n'); 
    fprintf('%s ', 'Orbit residuals: RMS(XYZ): ' ); 
    fprintf('%11.6f ', rms_orbc(1:3)); 
    fprintf('\n\n') 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Write orbit and partial derivatives to files in georb format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if write_data == 1 
    [OUT_data_foldername] = write_georb_data_products(orbit_config_struct, orbit_model_struct, orbit_matrices); 
else 
    OUT_data_foldername = ''; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
