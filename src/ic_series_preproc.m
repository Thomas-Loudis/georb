function [ic_data_struct] = ic_series_preproc(main_config_fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: ic_series_preproc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read configuration files and form the Initial Conditions array for all
%  sequential dates (initial epochs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - main_config_fname : Main Configuration file name
%
% Output arguments:
% - ic_data_struct    : Initial Conditions of all days in structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            21 January 2025
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

% Orbit time series :: No of sequential orbit arcs 
param_keyword = 'orbit_time_series_arcs';
[No_orbit_arcs_char] = read_param_file(main_config_fname,param_keyword);
No_orbit_arc_series = sscanf(No_orbit_arcs_char,'%d');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode: GEORB orbit mode (mission or objects) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode: georb_mode :: 'orbit_mission' : Orbit determination/propagation of satellite mission 
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

if  test_orbit_mission == 1
% Satellite missions cases :: GRACE missions
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission');

% GRACE or GRACE-FO mission' satellites
if test_mission_grace == 1
    satellite_1 = 'GRACE-A';
    satellite_2 = 'GRACE-B';
elseif test_mission_gracefo == 1
    satellite_1 = 'GRACE-C';
    satellite_2 = 'GRACE-D';
end

% Read IC configuration file to get the IC epochs list for the input orbiting mission/object
[ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_config_filename, satellite_1);
ic_data_satellite1 = ic_data_object1;

[ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_config_filename, satellite_2);
ic_data_satellite2 = ic_data_object1;

% Formulation of IC series matrix for N orbit arcs
[ic_sat1_series] = ic_series(ic_data_satellite1, No_orbit_arc_series);
[ic_sat2_series] = ic_series(ic_data_satellite2, No_orbit_arc_series);

ic_data_satellite1 = ic_sat1_series;
ic_data_satellite2 = ic_sat2_series;

% ic_data_struct(1).ic_data = ic_data_satellite1;
% ic_data_struct(2).ic_data = ic_data_satellite2;

% ic_series function outputs structure matrix matrixname.ic_data
ic_data_struct(1) = ic_data_satellite1;
ic_data_struct(2) = ic_data_satellite2;

else

% Mode: georb_mode :: 'orbit_objects' : Orbit determination of a single orbiting object (satellite, invidual orbiter)
% Read IC configuration file to get the IC epochs list for the input orbiting mission/object
[ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_config_filename, orbiting_object_name);

% Test the objects series option 'ALL_OBJECTS'
test_orbiting_object_name = strcmp(orbiting_object_name,'ALL_OBJECTS');
if test_orbiting_object_name == 1
    % ic_data_struct.ic_data = ic_data_objects;
    ic_data_objects_array = ic_data_objects;
else
    %ic_data_struct.ic_data = ic_data_object1;
    ic_data_objects_array = ic_data_object1;
end

% Fomrulation of IC series matrix for N orbit arcs
[ic_data_series] = ic_series(ic_data_objects_array, No_orbit_arc_series);
ic_data_struct = ic_data_series;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

