function [config_struct] = read_config2struct(config_filename) 

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
% Function:  read_param_cfg 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Read parameter' value of a configuration file *.IN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - config_filename : configuration file name 
% 
% Output arguments: 
% - config_struct   : Structure array of all configuration parameters  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                            27 October 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
fid = fopen(config_filename,'r'); 
i_struct = 0; 
while (~feof(fid)) 
    lineith = fgetl(fid);     
    param_name = sscanf(lineith,'%s %*'); 
    param_value = sscanf(lineith,'%*s %s %*'); 
    param_values_line = sscanf(lineith,'%*s %3000c'); 
    i_struct = i_struct + 1; 
    config_struct(i_struct).names  = param_name; 
    config_struct(i_struct).values = param_values_line; 
end 
fclose(fid); 
