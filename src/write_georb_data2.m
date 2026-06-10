function [georb_data_name] = write_georb_data2(orbit_config_fname, data_matrix, data_functional, reference_frame) 

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
% Function:  write_georb_data2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Write computed data to output files in GEORB format e.g. orbits, partial 
%  derivatives, observation residuals, external orbit comparison 
%  differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbit_config_fname  : GEORB configuration structure name  
% 
% Output arguments: 
% - georb_data_name     : GEORB written data file name 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                               9 July  2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 09/11/2022  Dr. Thomas Loudis Papanikolaou 
%             Minor update for introducing write_data_name function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Data functional name part1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% data_functional_part1 = 'orbit'; 
data_functional_part1 = data_functional; 
 
% Orbit 
data_functional_keyword = 'Orbit'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_part1 = 'orbit'; 
end 
 
% Partial Derivatives 
data_functional_keyword = 'Partial Derivatives'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_part1 = 'partials'; 
end 
 
% Orbit observation residuals 
data_functional_keyword = 'Orbit Observations residuals'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_part1 = 'obs_residuals'; 
end 
 
% Forces model: Acceleration Vector per epoch 
data_functional_keyword = 'Forces acceleration'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_part1 = 'aceleration'; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Data functional name part2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_functional_part2 = 'crf'; 
 
% ICRF 
reference_frame_keyword = 'ICRF'; 
test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
if test_reference_frame == 1 
data_functional_part2 = 'crf'; 
end 
 
% ITRF 
reference_frame_keyword = 'ITRF'; 
test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
if test_reference_frame == 1 
data_functional_part2 = 'trf'; 
end 
 
% ITRF 
reference_frame_keyword = 'Kepler'; 
test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
if test_reference_frame == 1 
data_functional_part2 = 'Kepler'; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Data functional file name 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_functional_filename = sprintf('%s%s%s%s','_',data_functional_part1,'_',data_functional_part2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% External Orbit Comparison differences  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_functional_keyword = 'External Orbit comparison'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
    % Orbital Frame 
    reference_frame_keyword = 'Orbital Frame'; 
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
    if test_reference_frame == 1 
        data_functional_filename = '_delta-RTN'; 
    end 
    % ICRF 
    reference_frame_keyword = 'ICRF'; 
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
    if test_reference_frame == 1 
        data_functional_filename = '_delta-CRF'; 
    end 
    % Kepler 
    reference_frame_keyword = 'Kepler'; 
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword); 
    if test_reference_frame == 1 
        data_functional_filename = '_delta-Kepler'; 
    end     
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inter-satellite observation residuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_functional_keyword = 'LRI range-rate residuals'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_filename = '_LRI_rangerate_residuals'; 
end 
 
data_functional_keyword = 'LRI range residuals'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_filename = '_LRI_range_residuals'; 
end 
 
data_functional_keyword = 'KBR range-rate residuals'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_filename = '_KBR_rangerate_residuals'; 
end 
 
data_functional_keyword = 'KBR range residuals'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_filename = '_KBR_range_residuals'; 
end 
 
data_functional_keyword = sscanf(data_functional,'%s%*'); 
test_data_functional = strcmp(data_functional_keyword,'Intersatellite'); 
if test_data_functional == 1 
data_functional_filename = '_intersatellite_ranging'; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Estimated Parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_functional_keyword = 'Parameters'; 
test_data_functional = strcmp(data_functional,data_functional_keyword); 
if test_data_functional == 1 
data_functional_filename = '_parameters'; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output Data file name to be written 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
mission_01 = 0; 
data_functional_keyword = sscanf(data_functional,'%3c%*'); 
test_data_functional = strcmp(data_functional_keyword,'LRI'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
data_functional_keyword = sscanf(data_functional,'%3c%*'); 
test_data_functional = strcmp(data_functional_keyword,'KBR'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
data_functional_keyword = sscanf(data_functional,'%8c%*'); 
test_data_functional = strcmp(data_functional_keyword,'intersat'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
 
data_functional_keyword = sscanf(data_functional,'%s%*'); 
test_data_functional = strcmp(data_functional_keyword,'Intersatellite'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
data_functional_keyword = sscanf(data_functional,'%s%*'); 
test_data_functional = strcmp(data_functional_keyword,'intersatellite'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
 
data_functional_keyword = sscanf(data_functional,'%3c%*'); 
test_data_functional = strcmp(data_functional_keyword,'NEQ'); 
if test_data_functional == 1 
mission_01 = 1; 
end 
% Output data file name conventions 
[georb_dataformat_name, georb_dataformat_suffix] = write_data_name(orbit_config_fname,mission_01); 
% Output Data file name 
georb_data_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GEORB version 
param_keyword = 'GEORB_version'; 
[src_version] = read_param_cfg(orbit_config_fname, param_keyword); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Write :: data matrix in georb format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[fid] = write_georb_format(georb_data_name, orbit_config_fname, data_matrix, data_functional, src_version, reference_frame); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
