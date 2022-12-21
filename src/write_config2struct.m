function [config_struct] = write_config2struct(georb_config_fname, orbit_model_fname, ic_data, src_version)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_orb_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write all input configurable variables to one configuration structure (file) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - georb_config_fname  : GEORB master configuration file name 
% - orbit_model_fname   : Orbit modelling configuration file name
% - ic_data             : Initial Conditions (IC) data line for the Orbiting-Object ID
% - src_version         : GEORB version
%
% Output arguments:
% - config_struct       : Strucutre array of all configuration varaibles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            18 August  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 16/11/2022, Dr. Thomas Loudis Papanikolaou
%             Modified to remove the satellite data configuration file and
%             move towards creating the data file names according to the
%             data format conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB general configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_struct = 0;

param_keyword = 'GEORB_version';
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = src_version;

param_keyword = 'georb_mode';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;

georb_mode = param_line;

param_keyword = 'orbiting_object_name';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
param_keyword_georb = 'orbiting_objects_mission';
config_struct(i_struct).names  = param_keyword_georb;
config_struct(i_struct).values = param_line;

orbiting_object_name = param_line;
mission_name = param_line;

param_keyword = 'orbit_time_series';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;

param_keyword = 'orbit_time_series_arcs';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;

param_keyword = 'ic_config_filename';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;

param_keyword = 'ic_state_vector_apriori';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;

% Orbit modelling configuration files
param_keyword = 'orb_config_filename';
[param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_line;
orbit_model_filename = param_value;

% % Satellite Data configuration file
% param_keyword = 'sat_data_filename';
% [param_value,param_line] = read_param_file(georb_config_fname,param_keyword);
% i_struct = i_struct + 1;
% config_struct(i_struct).names  = param_keyword;
% config_struct(i_struct).values = param_line;
% sat_data_filename = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode: georb_mode 
test_georb_mode = strcmp(georb_mode,'orbit_mission');
if test_georb_mode == 1
    mission_name = orbiting_object_name;
else
    satellite_id_name = sscanf(ic_data,'%s%*');
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-C');
    if test_satellite_mission == 1
       mission_name = 'GRACE_FO_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-D');
    if test_satellite_mission == 1
       mission_name = 'GRACE_FO_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-A');
    if test_satellite_mission == 1
       mission_name = 'GRACE_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-B');
    if test_satellite_mission == 1
       mission_name = 'GRACE_mission'; 
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Model configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config_file = -9;
orbit_config_fname = ' ';

if config_file == 1

fid = fopen(orbit_config_fname,'r');
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

end

if config_file == 0
   [config_struct_ORBIT] = read_config2struct(orbit_config_fname);    
   % Merge Structures
   [n,m] = size(config_struct_ORBIT);
   for i2_struct = 1 : n
        i_struct = i_struct + 1;
        config_struct(i_struct).names  = config_struct_ORBIT(i2_struct).names;
        config_struct(i_struct).values = config_struct_ORBIT(i2_struct).values;       
   end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Model configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Orbit Model file name to the input argument
orbit_model_filename = orbit_model_fname;

% Current orbit object/satellite name
orbiting_object_name = sscanf(ic_data,'%s %*');

param_keyword = 'orbiting_object_name';
param_value   = orbiting_object_name;
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'orbit_mode';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'Orbit_arc_length';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%param_keyword = orbiting_object_name;
%[IC_var_value, IC_dataline_c] = read_param_file(ic_data,param_keyword);
IC_dataline_c = sscanf(ic_data,'%*s %10000c %*');

% IC parameters line
param_keyword = 'IC_parameters';
param_value = sscanf(ic_data,'%*s %*s %*s %*s %*s %*s %*s %*s %*s %10000c %*');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Reference_frame';
param_value = sscanf(IC_dataline_c,'%s %*');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Time_scale';
param_value = sscanf(IC_dataline_c,'%*s %s %*');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Date_format';
param_value = sscanf(IC_dataline_c,'%*s %*s %s %*');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

Date_format = param_value;

param_keyword = 'Initial_Epoch';
param_date1 = sscanf(IC_dataline_c,'%*s %*s %*s %s %*');
param_date2 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %s %*');
param_date3 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %s %*');
param_date4 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %*s %s %*');
%param_value = sprintf('%s %s %s %s',param_date1, param_date2, param_date3, param_date4);
param_value = sprintf('%s %s %s %s',param_date4, param_date3, param_date2, param_date1);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


test = strcmp(Date_format,'MJD');
if test == 1
    MJD_day = sscanf(IC_dataline_c,'%*s %*s %*s %d %*');
    Sec_00h = sscanf(IC_dataline_c,'%*s %*s %*s %*s %f %*');
end

test = strcmp(Date_format,'calendar');
if test == 1
    year  = sscanf(IC_dataline_c,'%*s %*s %*s %d %*');
    month = sscanf(IC_dataline_c,'%*s %*s %*s %*s %d %*');
    day   = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %d %*');
    sec   = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %*s %f %*');
    [JD_ith, MJD_ith] = MJD_date(sec, day, month, year);
    MJD_day = fix(MJD_ith);
    Sec_00h = sec;
end

param_keyword = 'N_orbit_arcs';
param_value = sscanf(IC_dataline_c,'%*s %*s %*s  %*s %*s %*s %*s  %s %*');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'State_vector';
param_value = sscanf(IC_dataline_c,'%*s %*s %*s  %*s %*s %*s %*s  %*s  %10000c');
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation
param_keyword = 'EOP_filename';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'EOP_interpolation_points';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
% fprintf(fid,'%-27s %s ',param_keyword, param_value);
% fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'precession_nutation_model';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
% fprintf(fid,'%-27s %s ',param_keyword, param_value);
% fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration methods
param_keyword = 'Integration_method';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Stepsize';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'integrator_order_multistep';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Start_integrator';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'RKN_lamda_coefficient';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'RKN_interpolation_sigma';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational effects
param_keyword = 'Earth_Gravity_Field';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = '3rd_body_peturbations';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Tides_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'Relativity_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity field
param_keyword = 'Gravity_field_terms';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'gravity_model_fname';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'gravity_model_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'gravity_model_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


param_keyword = 'veq_gravity_model_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'veq_gravity_model_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planetary/Lunar ephemeris
param_keyword = 'planetary_ephemeris_DE_filename';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-35s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'planetary_ephemeris_DE_headername';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-35s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Tides
param_keyword = 'solid_earth_tides_1_non_freq';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'solid_earth_tides_2_freq';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides
param_keyword = 'ocean_tides';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'ocean_tides_model_fname';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'ocean_tides_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
%fprintf(fid,'%-27s %s ',param_keyword, param_value);
%fprintf(fid,'%s\n','');
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'ocean_tides_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'veq_ocean_tides_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'veq_ocean_tides_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pole Tide

% Solid Earth Pole Tide
param_keyword = 'solid_earth_pole_tide';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% Ocean Pole Tide
param_keyword = 'ocean_pole_tide';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric Tides
param_keyword = 'atm_tides';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'atm_tides_data_level';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'atm_tides_data_release';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'atm_tides_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% param_keyword = 'atm_tides_order';
% [param_value] = read_param_file(orbit_model_filename,param_keyword);
% i_struct = i_struct + 1;
% config_struct(i_struct).names  = param_keyword;
% config_struct(i_struct).values = param_value;
% 
% param_keyword = 'veq_atm_tides_degree';
% [param_value] = read_param_file(orbit_model_filename,param_keyword);
% i_struct = i_struct + 1;
% config_struct(i_struct).names  = param_keyword;
% config_struct(i_struct).values = param_value;
% 
% param_keyword = 'veq_atm_tides_order';
% [param_value] = read_param_file(orbit_model_filename,param_keyword);
% i_struct = i_struct + 1;
% config_struct(i_struct).names  = param_keyword;
% config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'aod_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AOD Data :: GRACE Level 1b data %AOD1B_2021-07-17_X_06.asc
% Create file name according to name format conventions based on data level name and date
param_keyword = 'aod_data_level';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
aod_data_level = param_value;

param_keyword = 'aod_data_release';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
aod_data_release = param_value;

%datalevel_name  = 'AOD1B';
%datarelease_ver = '06';

datalevel_name  = aod_data_level;
datarelease_ver = aod_data_release;
% Data file name based on data naming conventions
[AOD_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);

% Write file name to the orbit configuration file 
param_keyword = 'AOD_data_filename';
param_value = AOD_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_keyword = 'aod_effect_data_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'AOD_degree_max';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic Effects
param_keyword = 'Schwarzschild_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'Lense_Thirring_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'geodesic_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PPN_beta_parameter';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PPN_gama_parameter';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'C_speed_of_light';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Gravitational effects
param_keyword = 'non_gravitational_forces';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-A   calendar 2009 11 17   GNV1B_2009-11-17_A_01.asc    GNV1B_2009-11-17_A_01.asc      ACC1B_2009-11-17_A_01.asc      SCA1B_2009-11-17_A_01.asc 
% GRACE-C   calendar 2021 07 12    GRACEFO-1_kinematicOrbit_2021-07-12.txt	GNV1B_2021-07-12_C_04.txt	ACT1B_2021-07-12_C_04.txt   SCA1B_2021-07-12_C_04.txt   KBR1B_2021-07-12_Y_04.txt   LRI1B_2021-07-12_Y_04.txt
% id1_satellite = orbiting_object_name;
% id2_MJD = MJD_day;
% % Satellite Data line for the current MJD day 
% [satdata_line_mjd] = read_satdata_cfg(sat_data_filename,id1_satellite,id2_MJD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'acc_data';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'acc_cal_paramestim';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: Accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%param_keyword = 'accelerometer_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s %s %*');

param_keyword = 'acc_data_level';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
acc_data_level = param_value;

param_keyword = 'acc_data_release';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
acc_data_release = param_value;

%datalevel_name  = 'ACT1B';
%datarelease_ver = '04';
datalevel_name  = acc_data_level;
datarelease_ver = acc_data_release;

% Accelerometer dat file name according to name conventions and input date
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);

param_keyword = 'accelerometer_data';
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: Star Camera data
%param_keyword = 'star_camera_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s %s %*');

param_keyword = 'sca_data_level';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
sca_data_level = param_value;

param_keyword = 'sca_data_release';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
sca_data_release = param_value;

%datalevel_name  = 'SCA1B';
%datarelease_ver = '04';
datalevel_name  = sca_data_level;
datarelease_ver = sca_data_release;

% Star Camera data file name according to name conventions and input date
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);

param_keyword = 'star_camera_data';
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer data processing parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'accelerometer_interp_no';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'star_camera_interp_no';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer Calibration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer Calibration modelling
param_keyword = 'acc_cal_bias';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'acc_cal_bias_drift_1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'acc_cal_bias_drift_2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'acc_cal_scale';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'empirical_forces';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'empirical_frame';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% Bias terms
param_keyword = 'empirical_bias_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'empirical_bias_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'empirical_bias_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% CPR terms
param_keyword = 'cpr_C_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_S_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_C_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_S_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_C_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_S_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'cpr_freq_number';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------------------------------------------------------------
% Empirical Acclerations / Pulses (Piecewise accelerations or Velocity changes)
% ---------------------------------------------------------------------------
param_keyword = 'PULSES_estim_yn';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
  
param_keyword = 'stoch_param_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_epochs_number';
% [param_value] = read_param_file(orbit_model_filename,param_keyword);
param_value = '000';
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_interval';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'stoch_time_interval';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
    
param_keyword = 'PULSES_offset';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_frame';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_axis_1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_axis_2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'PULSES_axis_3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
% ---------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'pseudo_obs_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
pseudo_obs_type = param_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: Pseudo-Observation data
% GRACE-A   calendar 2009 11 17   GNV1B_2009-11-17_A_01.asc    GNV1B_2009-11-17_A_01.asc   ACC1B_2009-11-17_A_01.asc      SCA1B_2009-11-17_A_01.asc 
param_keyword = 'pseudo_obs_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s %s %*');
datalevel_name  = 'kinematicOrbit';
obs_keyword = 'itsg_kin';
test_obs_type = strcmp(pseudo_obs_type, obs_keyword);
if test_obs_type == 1
datalevel_name  = 'kinematicOrbit';
end
obs_keyword = 'gnv1b';
test_obs_type = strcmp(pseudo_obs_type, obs_keyword);
if test_obs_type == 1
datalevel_name  = 'GNV1B';
end
%datarelease_ver = '04';
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_keyword = 'cov_pseudo_obs_data';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'obs_outliers_yn';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'pseudo_obs_sigma';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Satellite Ranging Observations :: GRACE & GRACE-Follow On missions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: KBR data
param_keyword = 'KBR_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s%*s %s %*');
datalevel_name  = 'KBR1B';
datarelease_ver = '04';
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% Satellite Data :: LRI data
param_keyword = 'LRI_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s %s %*');
datalevel_name  = 'LRI1B';
datarelease_ver = '04';
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimator algorithm
param_keyword = 'estimator_iterations';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'external_orbit_comp';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

param_keyword = 'external_orbit_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
external_orbit_type = param_value;

% External Orbit data format ('kinematic' or 'dynamic'):
% Options:
% - 'gnv1b'         : GRACE mission GNV1b orbit data format 
% - 'georb'         : Internal orbit data format of GEORB software
% - 'kepler'        : Kepler orbit :: "Status: Not supported by current version"
% external_orbit_type            gnv1b 

% Satellite Data file :: External orbit data
param_keyword = 'external_orbit_data';
%param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s %s %*');
datalevel_name  = 'GNV1B';
extorb_keyword = 'gnv1b';
test_obs_type = strcmp(external_orbit_type, extorb_keyword);
if test_obs_type == 1
datalevel_name  = 'GNV1B';
end
%datarelease_ver = '04';
[data_filename_conv] = data_name_conv(mission_name, orbiting_object_name, datalevel_name, datarelease_ver, MJD_day);
param_value = data_filename_conv;
i_struct = i_struct + 1;
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
