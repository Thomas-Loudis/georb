function [fid] = write_orb_config(out_config_filename, orbiting_object_name, ic_data, orbit_model_filename, sat_data_filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_orb_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write all input configurable variables to one orbit configuration file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - out_filename  : Orbit output file name
% - cfg_filename  : Configuration file name
%
% Output arguments:
% -     : 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            18 August  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Open file for writing
fid = fopen(out_config_filename,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = '---------------------------------------------------------------------------';
fprintf(fid,'%s \n',param_keyword);

param_keyword = 'Orbiting object configutation file';
fprintf(fid,'%s \n',param_keyword);

param_keyword = '---------------------------------------------------------------------------';
fprintf(fid,'%s \n',param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write configuration variables

param_keyword = 'orbiting_object_name';
param_value   = orbiting_object_name;
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'orbit_mode';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Orbit_arc_length';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-C     ICRF   GPS   calendar   2019 07 18   0.0     1    3317114.29453071626  252279.97141062448 -6033407.00916356966 -6639.98945932193874 -227.41515344014149 -3673.76634420864957
%param_keyword = orbiting_object_name;
%[IC_var_value, IC_dataline_c] = read_param_file(ic_data,param_keyword);
IC_dataline_c = sscanf(ic_data,'%*s %10000c %*');

% IC parameters line
param_keyword = 'IC_parameters';
param_value = sscanf(ic_data,'%*s %*s %*s %*s %*s %*s %*s %*s %*s %10000c %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


param_keyword = 'Reference_frame';
param_value = sscanf(IC_dataline_c,'%s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Time_scale';
param_value = sscanf(IC_dataline_c,'%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Date_format';
param_value = sscanf(IC_dataline_c,'%*s %*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
Date_format = param_value;

param_keyword = 'Initial_Epoch';
param_date1 = sscanf(IC_dataline_c,'%*s %*s %*s %s %*');
param_date2 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %s %*');
param_date3 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %s %*');
param_date4 = sscanf(IC_dataline_c,'%*s %*s %*s %*s %*s %*s %s %*');
%param_value = sprintf('%s %s %s %s',param_date1, param_date2, param_date3, param_date4);
param_value = sprintf('%s %s %s %s',param_date4, param_date3, param_date2, param_date1);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

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
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'State_vector';
param_value = sscanf(IC_dataline_c,'%*s %*s %*s  %*s %*s %*s %*s  %*s  %10000c');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth Orientation
param_keyword = 'EOP_filename';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'EOP_interpolation_points';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'precession_nutation_model';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Numerical Integration methods
param_keyword = 'Integration_method';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Stepsize';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'integrator_order_multistep';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Start_integrator';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'RKN_lamda_coefficient';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'RKN_interpolation_sigma';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Gravitational effects
param_keyword = 'Earth_Gravity_Field';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = '3rd_body_peturbations';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Tides_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Relativity_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Gravity field
param_keyword = 'Gravity_field_terms';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'gravity_model_fname';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'gravity_model_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'gravity_model_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'veq_gravity_model_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'veq_gravity_model_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Planetary/Lunar ephemeris
param_keyword = 'planetary_ephemeris_DE_filename';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-35s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'planetary_ephemeris_DE_headername';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-35s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Solid Earth Tides
param_keyword = 'solid_earth_tides_1_non_freq';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'solid_earth_tides_2_freq';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Ocean Tides
param_keyword = 'ocean_tides';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'ocean_tides_model_fname';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'ocean_tides_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'ocean_tides_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'veq_ocean_tides_degree';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'veq_ocean_tides_order';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AOD Data :: GRACE Level 1b data
% Create file name according to name format conventions based on data level name and date
%AOD1B_2021-07-17_X_06.asc
%MJD_day = 59412

data_name_level = 'AOD1B';
release_ver     = '06';
format_ext      = '.asc';
[sec,day_no,month_no,year] = MJD_inv(MJD_day);

if day_no < 10
    fname_day = sprintf('%s%d','0',day_no);
else
    fname_day = sprintf('%d',day_no);
end    

if month_no < 10
    fname_month = sprintf('%s%d','0',month_no);
else
    fname_month = sprintf('%d',month_no);
end

% AOD file name considering format name conventions
AOD_filename_conv = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_name_level,'_',year,'-',fname_month,'-',fname_day,'_','X','_',release_ver,format_ext);

% Write file name to the orbit configuration file 
param_keyword = 'AOD_data_filename';
param_value = AOD_filename_conv;
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'aod_effects';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'aod_effect_data_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'AOD_degree_max';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Relativistic Effects
param_keyword = 'Schwarzschild_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'Lense_Thirring_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'geodesic_effect';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PPN_beta_parameter';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PPN_gama_parameter';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'C_speed_of_light';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


% Non-Gravitational effects
param_keyword = 'non_gravitational_forces';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-A   calendar 2009 11 17   GNV1B_2009-11-17_A_01.asc    GNV1B_2009-11-17_A_01.asc      ACC1B_2009-11-17_A_01.asc      SCA1B_2009-11-17_A_01.asc 
id1_satellite = orbiting_object_name;
id2_MJD = MJD_day;
% Satellite Data line for the current MJD day 
[satdata_line_mjd] = read_satdata_cfg(sat_data_filename,id1_satellite,id2_MJD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: Accelerometer data
param_keyword = 'accelerometer_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% Satellite Data :: Star Camera data
param_keyword = 'star_camera_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accelerometer data processing parameters
param_keyword = 'accelerometer_interp_no';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'star_camera_interp_no';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% Accelerometer Calibration modelling
param_keyword = 'acc_cal_bias';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'acc_cal_bias_drift_1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'acc_cal_bias_drift_2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'acc_cal_scale';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'empirical_forces';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'empirical_frame';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% Bias terms
param_keyword = 'empirical_bias_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'empirical_bias_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'empirical_bias_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% CPR terms
param_keyword = 'cpr_C_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_S_axis1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_C_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_S_axis2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_C_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_S_axis3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cpr_freq_number';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------------------------------------------------------------
% Stochastic Pulses (velocity or acceleration changes)
% ---------------------------------------------------------------------------
param_keyword = 'PULSES_estim_yn';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
  
param_keyword = 'stoch_param_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_epochs_number';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_interval';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'stoch_time_interval';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
    
param_keyword = 'PULSES_offset';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_frame';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_axis_1';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_axis_2';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'PULSES_axis_3';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
% ---------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: Pseudo-Observation data
% GRACE-A   calendar 2009 11 17   GNV1B_2009-11-17_A_01.asc    GNV1B_2009-11-17_A_01.asc   ACC1B_2009-11-17_A_01.asc      SCA1B_2009-11-17_A_01.asc 
param_keyword = 'pseudo_obs_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'pseudo_obs_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'cov_pseudo_obs_data';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'pseudo_obs_sigma';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Satellite Ranging Observations :: GRACE & GRACE-Follow On missions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Data :: KBR data
param_keyword = 'KBR_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% Satellite Data :: LRI data
param_keyword = 'LRI_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimator algorithm
param_keyword = 'estimator_iterations';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'external_orbit_comp';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

% Satellite Data file :: External orbit data
param_keyword = 'external_orbit_data';
param_value = sscanf(satdata_line_mjd,'%*s%*s%*s%*s%*s%*s %s %*');
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');

param_keyword = 'external_orbit_type';
[param_value] = read_param_file(orbit_model_filename,param_keyword);
fprintf(fid,'%-27s %s ',param_keyword, param_value);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fid);
