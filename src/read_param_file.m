function [param_value, param_line] = read_param_file(config_name,param_keyword) 

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
% Function:  read_param_file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Read parameter' value of a configuration file *.IN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - config_name     : prm.in input file name (configuration file) 
% - param_keyword  : Parameter keyword  
% 
% Output arguments: 
% - param_value    : Parameter value from the configuration file 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Papanikolaou                                       27 May 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 27/10/2022  Dr. Thomas Loudis Papanikolaou 
%             Options for the orbit configuration format based on structure array or file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
config_mode = 1; 
 
% Orbit configuration based on file 
if config_mode == 1 
param_value = ''; 
param_line = ''; 
fid = fopen(config_name,'r'); 
while (~feof(fid)) 
    lineith = fgetl(fid);     
    line_keyword = sscanf(lineith,'%s %*'); 
    test_keyword = strcmp(line_keyword,param_keyword); 
    if test_keyword == 1  
        param_value = sscanf(lineith,'%*s %s %*'); 
        param_line = sscanf(lineith,'%*s %3000c');  
    end         
end 
fclose(fid); 
end 
 
 
% Orbit configuration based on structure array 
if config_mode == 2 
[k m] = size(config_name); 
for i = 1 : m 
    param_name = config_name(1,i).names; 
    test_keyword = strcmp(param_name,param_keyword); 
    if test_keyword == 1  
        param_line = config_name(1,i).values; 
        param_value = sscanf(param_line,'%s %*'); 
    end             
end    
end 
