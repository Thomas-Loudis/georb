function [georb_dataformat_name, georb_dataformat_suffix] = write_data_name(orbit_config_fname, mission_01)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_data_name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write orbits and partial derivatives to data files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_config_fname  : GEORB orbit configuration structure name
% - mission_01          : Satellite mission flag
%
% Output arguments:
% -     : 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            9 November 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode: georb_mode :: 'orbit_mission' : Orbit determination of a single orbiting object (satellite, invidual orbiter)
param_keyword = 'georb_mode';
[georb_mode] = read_param_cfg(orbit_config_fname,param_keyword);
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_objects_mission';
[orbiting_object_mission_name] = read_param_cfg(orbit_config_fname,param_keyword);

% Satellite missions cases ::
test_grace   = strcmp(orbiting_object_mission_name,'GRACE_mission');
test_gracefo = strcmp(orbiting_object_mission_name,'GRACE_FO_mission');

% Satellite/Object name
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname,param_keyword);

% Read cofig file for MJDo
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data file name conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data format file name
if mission_01 == 1
% Case of satellite mission
    if test_grace == 1
        georb_dataformat_name = sprintf('%s%s%d', 'grace','_',fix(IC_MJDo));
    elseif test_gracefo == 1
        georb_dataformat_name = sprintf('%s%s%d', 'grace-fo','_',fix(IC_MJDo));
    end

else
    % Case of individual orbiting object/satellite
    georb_dataformat_name = sprintf('%s%s%d', orbiting_object_name,'_',fix(IC_MJDo));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data format suffix
georb_dataformat_suffix = sprintf('%s','.orb');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
