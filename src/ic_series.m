function [ic_time_series] = ic_series(ic_objects, No_orbit_arc_series) 

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
% Function: ic_series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  orbit_mission_grace is the main function for calling the source code 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - config_filename     : Configuration file name 
% 
% Output arguments: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                            25 August  2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ic_time_series (1,:) = ic_objects(1,:); 
[ic_n, ic_m] = size(ic_objects); 
ic_k = 0; 
N_orbit_arcs = No_orbit_arc_series; 
for ic_i = 1 : ic_n 
    % Current IC data line 
    ic_data_object_i = ic_objects(ic_i,:); 
    % Number of orbit arcs 
    % N_orbit_arcs = sscanf(ic_data_object_i,'%*s%*s%*s%*s%*s%*s%*s%*s %d %*') 
    % Create new IC data lines for the times series of N orbit arcs  
    for i_arc = 1 : N_orbit_arcs + 1 
        arc_length = 24; 
        [ic_data_object_Narcs] = ic_cfg_series(ic_data_object_i, i_arc, arc_length);        
        ic_k = ic_k + 1; 
        % ic_time_series (ic_k , :) = ic_data_object_Narcs; 
        Nc = length(ic_data_object_Narcs); 
        % ic_time_series.ic_data(ic_k , 1:Nc) = ic_data_object_Narcs; 
        ic_time_series(ic_k).ic_data = ic_data_object_Narcs;         
    end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
