function [ic_time_series] = ic_series(ic_objects)

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
ic_time_series (1,:) = ic_objects(1,:);
[ic_n ic_m] = size(ic_objects);
ic_k = 1;
for ic_i = 2 : ic_n
    % Current IC data line
    ic_data_object_i = ic_objects(ic_i,:);
    % Number of orbit arcs
    N_orbit_arcs = sscanf(ic_data_object_i,'%*s%*s%*s%*s%*s%*s%*s%*s %d %*');
    % Create new IC data lines for the times series of N orbit arcs 
    for i_arc = 1 : N_orbit_arcs        
        [ic_data_object_Narcs] = ic_cfg_series(ic_data_object_i, i_arc, arc_length);        
        ic_k = ic_k + 1;
        ic_time_series (ic_k , :) = ic_data_object_Narcs;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
