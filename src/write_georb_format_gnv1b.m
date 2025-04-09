function [fid] = write_georb_format_gnv1b(out_filename, cfg_filename, data_matrix, data_functional, georb_version,reference_frame)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_georb_format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write coputations to output file according to the GEORB data format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - out_filename  : Orbit output file name
% - cfg_filename  : Configuration file name
%
% Output arguments:
% -     : 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               9 July  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Open file for writing
fid = fopen(out_filename,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'GEORB format file';
fprintf(fid,'%-33s',param_keyword);
fprintf(fid,'%s\n','');

param_keyword = 'Data functional';
fprintf(fid,'%-33s %s %s',param_keyword,': ',data_functional);
fprintf(fid,'%s\n','');

param_keyword = 'Software_version';
fprintf(fid,'%-33s %s %s %s',param_keyword, ': ', 'GEORB', georb_version);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'Orbiting Object/Satellite';
fprintf(fid,'%-33s %s %s ',param_keyword, ': ', orbiting_object_name);
fprintf(fid,'%s\n','');

param_keyword = 'Reference Frame';
fprintf(fid,'%-33s %s %s ',param_keyword, ': ', reference_frame);
fprintf(fid,'%s\n','');

param_keyword = 'Time scale';
fprintf(fid,'%-33s %s %s ',param_keyword, ': ', 'GPS Time');
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date
mjd_to = data_matrix(1,1);
[UT,day,month,year] = MJD_inv(mjd_to);
param_keyword = 'Date';
fprintf(fid,'%-33s %s %d %d %d ',param_keyword, ': ',day,month,year);
fprintf(fid,'%s\n','');

[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_filename);
orbit_arc_hours = orbit_arc_length / 3600;
param_keyword = 'Orbit arc length (hours)';
fprintf(fid,'%-33s %s %f ',param_keyword, ': ', orbit_arc_hours);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration methods    
[integr_method, integr_stepsize, integr_order, integr_start, integr_RKN_lamda, integr_RKN_sigma] = prm_integrator(cfg_filename);

param_keyword = 'Numerical Integration';
fprintf(fid,'%-33s %s %s %d %s ',param_keyword, ': ', integr_method, integr_order, 'order');
fprintf(fid,'%s\n','');
% Orbit integrator step (sec)
param_keyword = 'Integrator step (sec)';
fprintf(fid,'%-33s %s %d ',param_keyword, ': ', integr_stepsize);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation                                                       
param_keyword = 'EOP_filename';
[eop_name] = read_param_cfg(cfg_filename,param_keyword);

param_keyword = 'Earth Orientation EOP';
fprintf(fid,'%-33s %s %s ',param_keyword, ': ', eop_name);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precession-Nutation model                                                       
param_keyword = 'precession_nutation_model';
[precession_nutation_model] = read_param_cfg(cfg_filename,param_keyword);

param_keyword = 'Precession Nutation model';
fprintf(fid,'%-33s %s %s ',param_keyword, ': ', precession_nutation_model);
fprintf(fid,'%s\n','');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Dynamics
[gravitational_effects_01, non_gravitational_forces_01, tidal_effects_01, Gravity_field_terms] = prm_dynamics(cfg_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model d/o : GOCO06s  100x100
param_keyword = 'gravity_model_fname';
[gfm_name] = read_param_cfg(cfg_filename,param_keyword);

param_keyword = 'gravity_model_degree';
[gravity_model_degree] = read_param_cfg(cfg_filename,param_keyword);
N_max_grav = sscanf(gravity_model_degree,'%d %*');

param_keyword = 'gravity_model_order';
[gravity_model_order] = read_param_cfg(cfg_filename,param_keyword);
M_max_grav = sscanf(gravity_model_order,'%d %*');

param_keyword = 'Gravity Field model d/o';
fprintf(fid,'%-33s %s %s %d%s%d',param_keyword, ': ', gfm_name, N_max_grav,'x',M_max_grav);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planetary_ephemeris     : -
if gravitational_effects_01(2,1) == 0
    %param_keyword = 'Third Body perturbations';
    param_keyword =  'Planetary/Lunar Ephemeris';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');
    fprintf(fid,'%s\n','');
else
    param_keyword = 'planetary_ephemeris_DE_filename';
    [param_value] = read_param_cfg(cfg_filename,param_keyword);
    Nchar = length(param_value);
    DE_no = param_value(Nchar-2 : Nchar);
    DE_name = sprintf('%s%s','DE',DE_no);

    param_keyword =  'Planetary/Lunar Ephemeris';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ',DE_name);
    fprintf(fid,'%s\n','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Tides       : - Frequency dependent Terms  
param_keyword = 'solid_earth_tides_2_freq';
[param_value] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'Solid Earth Tides';
fprintf(fid,'%-33s %s %s', param_keyword,': ',param_value);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides model d/o   100x100 
param_keyword = 'ocean_tides';
[param_value] = read_param_cfg(cfg_filename,param_keyword);

if param_value == 'y' 
    param_keyword = 'ocean_tides_model_fname';
    [octides_model_name] = read_param_cfg(cfg_filename,param_keyword);

    param_keyword = 'ocean_tides_degree';
    [ocean_tides_degree] = read_param_cfg(cfg_filename,param_keyword);
    N_max_otides = sscanf(ocean_tides_degree,'%d ');

    param_keyword = 'ocean_tides_order';
    [ocean_tides_order] = read_param_cfg(cfg_filename,param_keyword);
    M_max_otides = sscanf(ocean_tides_order,'%d ');

    param_keyword = 'Ocean Tides model d/o';
    fprintf(fid,'%-33s %s %s %d%s%d',param_keyword, ': ', octides_model_name, N_max_otides,'x',M_max_otides);
    %fprintf(fid,'%-33s %s %s %d%s%d',param_keyword, ': ', octides_model_name, ocean_Nmax_glb,'x',ocean_Mmax_glb); 
    fprintf(fid,'%s\n','');
else   
    param_keyword = 'Ocean Tides model d/o';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');    
    fprintf(fid,'%s\n','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Pole Tide
param_keyword = 'solid_earth_pole_tide';
[effect_yn] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'Solid Earth Pole Tide';
test = strcmp(effect_yn,'y');
if test == 0
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');    
    fprintf(fid,'%s\n','');    
elseif test == 1
    fprintf(fid,'%-33s %s %s', param_keyword,': ','IERS Conventions 2010');
    fprintf(fid,'%s\n','');
end

% Ocean Pole Tide
param_keyword = 'ocean_pole_tide';
[effect_yn] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'Ocean Pole Tide';
test = strcmp(effect_yn,'y');
if test == 0
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');    
    fprintf(fid,'%s\n','');    
elseif test == 1
    fprintf(fid,'%-33s %s %s', param_keyword,': ','IERS Conventions 2010');
    fprintf(fid,'%s\n','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric Tides 
param_keyword = 'atm_tides';
[atm_tides_effects_yn] = read_param_cfg(cfg_filename,param_keyword);
test = strcmp(atm_tides_effects_yn,'y');
if test == 0
    param_keyword = 'Atmospheric Tides'; 
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');    
    fprintf(fid,'%s\n','');    
elseif test == 1
    param_keyword = 'atm_tides_data_level';
    [atm_tides_data_level] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'atm_tides_data_release';
    [atm_tides_data_release] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'atm_tides_degree';
    [atm_tides_degree] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'Atmospheric Tides'; 
    fprintf(fid,'%-33s %s %s%s %s %s', param_keyword,': ',atm_tides_data_level,atm_tides_data_release, 'Degree', atm_tides_degree);
    fprintf(fid,'%s\n','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmosphere and Ocean De-Aliasing effects
param_keyword = 'aod_effects';
[aod_effects_yn] = read_param_cfg(cfg_filename,param_keyword);
test = strcmp(aod_effects_yn,'y');
if test == 0
    param_keyword = 'Atmosphere and Ocean De-Aliasing'; 
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');    
    fprintf(fid,'%s\n','');    
elseif test == 1
    param_keyword = 'AOD_data_filename';
    [AOD_data_filename] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'Atmosphere and Ocean De-Aliasing'; 
    fprintf(fid,'%-33s %s %s', param_keyword,': ',AOD_data_filename);
    fprintf(fid,'%s\n','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic effects
param_keyword = 'Relativistic effects';

if gravitational_effects_01(4,1) == 0
fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');
fprintf(fid,'%s\n','');
else
fprintf(fid,'%-33s %s %s', param_keyword, ': ', 'IERS Conventions 2010');
fprintf(fid,'%s\n','');        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Gravitational Forces

% Accelerometry
if non_gravitational_forces_01 == 0
    param_keyword = 'Accelerometer data';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', 'n');
    fprintf(fid,'%s\n','');
else
    % Number of Calibration parameters estimated
    % Bias and Bias' Drfit terms
    Nparam_ACC_CAL = 0;
    param_keyword = 'acc_cal_bias';
    [cal_param] = read_param_cfg(cfg_filename,param_keyword);
    if cal_param == 'y'
    Nparam_ACC_CAL = Nparam_ACC_CAL + 3;
    end
    param_keyword = 'acc_cal_bias_drift_1';
    [cal_param] = read_param_cfg(cfg_filename,param_keyword);
    if cal_param == 'y'
    Nparam_ACC_CAL = Nparam_ACC_CAL + 3;
    end
    param_keyword = 'acc_cal_bias_drift_2';
    [cal_param] = read_param_cfg(cfg_filename,param_keyword);
    if cal_param == 'y'
    Nparam_ACC_CAL = Nparam_ACC_CAL + 3;
    end
    % Scale factors
    param_keyword = 'acc_cal_scale';
    [cal_param] = read_param_cfg(cfg_filename,param_keyword);
    test = strcmp(cal_param,'diagonal');
    if test == 1
        Nparam_ACC_CAL = Nparam_ACC_CAL + 3;
    end
    test = strcmp(cal_param,'semi-full');
    if test == 1
        Nparam_ACC_CAL = Nparam_ACC_CAL + 6;
    end
    test = strcmp(cal_param,'full');
    if test == 1
        Nparam_ACC_CAL = Nparam_ACC_CAL + 9;
    end
    
    % Accelerometer data file
    param_keyword = 'accelerometer_data';
    [accelerometer_data_file] = read_param_cfg(cfg_filename,param_keyword);
    
    % Star Camera data file
    param_keyword = 'star_camera_data';
    [star_camera_data_file] = read_param_cfg(cfg_filename,param_keyword);
    
    % Write to file
    param_keyword = 'Accelerometer Data';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', accelerometer_data_file);
    fprintf(fid,'%s\n','');        

    param_keyword = 'Star Camera Data';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', star_camera_data_file);
    fprintf(fid,'%s\n','');        

    param_keyword = 'Accelerometer Calibration';
    fprintf(fid,'%-33s %s %s %d %s', param_keyword, ': ', 'Calibration model', Nparam_ACC_CAL,'parameters');
    fprintf(fid,'%s\n','');        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces
param_keyword = 'empirical_forces';
[empirical_forces_yn] = read_param_cfg(cfg_filename,param_keyword);

if empirical_forces_yn == 'n'
param_keyword = 'Empirical Accelerations (1)';
fprintf(fid,'%-33s %s %s', param_keyword, ': ', 'n');
fprintf(fid,'%s\n','');    
elseif empirical_forces_yn == 'y'

param_keyword = 'empirical_bias_axis1';
[empirical_bias_axis1] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'empirical_bias_axis2';
[empirical_bias_axis2] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'empirical_bias_axis3';
[empirical_bias_axis3] = read_param_cfg(cfg_filename,param_keyword);
Np_bias = 0;
if empirical_bias_axis1 == 'y'
Np_bias = Np_bias + 1;
end
if empirical_bias_axis2 == 'y'
Np_bias = Np_bias + 1;
end
if empirical_bias_axis3 == 'y'
Np_bias = Np_bias + 1;
end

Np_cpr = 0;
param_keyword = 'cpr_C_axis1';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end
param_keyword = 'cpr_S_axis1';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end
param_keyword = 'cpr_C_axis2';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end
param_keyword = 'cpr_S_axis2';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end
param_keyword = 'cpr_C_axis3';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end
param_keyword = 'cpr_S_axis3';
[cpr_term] = read_param_cfg(cfg_filename,param_keyword);
if cpr_term == 'y'
Np_cpr = Np_cpr + 1;
end

param_keyword = 'empirical_frame';
[empirical_frame] = read_param_cfg(cfg_filename,param_keyword);

param_keyword = 'Empirical Accelerations (1)';
fprintf(fid,'%-33s %s %s %s %s %d %s %d', param_keyword, ': ','Frame:', empirical_frame, 'Bias:', Np_bias ,'Cycle-per-revolution:', Np_cpr);
fprintf(fid,'%s\n','');        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces :: Pulses or Piecewise accelerations
param_keyword = 'PULSES_estim_yn';
[PULSES_estim_yn] = read_param_cfg(cfg_filename,param_keyword);

test = strcmp(PULSES_estim_yn,'y');
if test == 0
    param_keyword = 'Empirical Accelerations (2)';
    fprintf(fid,'%-33s %s %s', param_keyword, ': ', '-');
    fprintf(fid,'%s\n','');
    
elseif test == 1    
    param_keyword = 'stoch_param_type';
    [stoch_param_type] = read_param_cfg(cfg_filename,param_keyword);
    test_stoch = strcmp(stoch_param_type,'stoch_pulses');
    if test_stoch == 1
        stoch_param_type = 'Pulses';
    end
    test_stoch = strcmp(stoch_param_type,'stoch_accel_constant');
    if test_stoch == 1
        stoch_param_type = 'Piecewise constant accelerations';
    end
    
    param_keyword = 'PULSES_frame';
    [empirical_frame] = read_param_cfg(cfg_filename,param_keyword);

    param_keyword = 'PULSES_interval';
    [PULSES_interval] = read_param_cfg(cfg_filename,param_keyword);
    Step_interval_hours = sscanf(PULSES_interval, '%d%*') / 3600;
    
    param_keyword = 'stoch_time_interval';
    [stoch_time_interval] = read_param_cfg(cfg_filename,param_keyword);
    duration_min = sscanf(stoch_time_interval, '%d%*') / 60;
    
    param_keyword = 'PULSES_axis_1';
    [PULSES_axis_1] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'PULSES_axis_2';
    [PULSES_axis_2] = read_param_cfg(cfg_filename,param_keyword);
    param_keyword = 'PULSES_axis_3';
    [PULSES_axis_3] = read_param_cfg(cfg_filename,param_keyword);
    
    param_keyword = 'Empirical Accelerations (2)';
    fprintf(fid,'%-33s %s %s %s %s%s%s%s %s %d %s %d', param_keyword, ': ','Frame:', empirical_frame, 'Directions Vector: ', PULSES_axis_1,PULSES_axis_2,PULSES_axis_3, 'Step interval (hours):', Step_interval_hours ,'Duration (min):', duration_min);
    fprintf(fid,'%s\n','');        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-Observations
if 1 <0
param_keyword = 'pseudo_obs_type';
[pseudo_obs_type] = read_param_cfg(cfg_filename,param_keyword);
param_keyword = 'Observations';
fprintf(fid,'%-33s %s %s %s', param_keyword, ': ', 'Pseudo-Observations orbit data: ' ,pseudo_obs_type);
fprintf(fid,'%s\n','');

% Pseudo-observations orbit data file 
%pseudo_obs_data      =     GRACEFO-2_kinematicOrbit_2019-07-18.txt 

param_keyword = 'pseudo_obs_data';
[pseudo_obs_data] = read_param_cfg(cfg_filename,param_keyword) 
param_keyword = 'Pseudo-Observations orbit data';
fprintf(fid,'%-33s %s %', param_keyword, ': ',pseudo_obs_data);
fprintf(fid,'%s\n','');
end

param_keyword = 'pseudo_obs_type';
[pseudo_obs_type] = read_param_cfg(cfg_filename,param_keyword);

param_keyword = 'pseudo_obs_data';
[pseudo_obs_data] = read_param_cfg(cfg_filename,param_keyword); 

param_keyword = 'Observations';
fprintf(fid,'%-33s %s %s %s %s', param_keyword, ': ', pseudo_obs_type, ' orbit data: ', pseudo_obs_data);
fprintf(fid,'%s\n','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data lines format             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format_line = 'Modified Julian Day number | Seconds since 00h | X(m) Y(m) Z(m) | Vx(m/sec) Vy(m/sec) Vz(m/sec)' ;

% Orbit Format
data_product = 'Orbit';
test_data_functional = strcmp(data_functional,data_product);
if test_data_functional == 1
format_line = 'Modified Julian Day number | Seconds since 00h | X(m) Y(m) Z(m) | Vx(m/sec) Vy(m/sec) Vz(m/sec)' ;
end

% Pseudo-Observation residuals
data_product  = 'Orbit Observations residuals';
test_data_functional = strcmp(data_functional,data_product);
if test_data_functional == 1
format_line = 'Modified Julian Day number | Seconds since 00h | X(m) Y(m) Z(m)' ;
end

% Partials
data_product  = 'Partial Derivatives';
test_data_functional = strcmp(data_functional,data_product);
if test_data_functional == 1
format_line = 'Modified Julian Day number | Seconds since 00h | dZ(i)/dXo dZ(i)/dYo dZ(i)/dZo dZ(i)/dVxo dZ(i)/dVyo dZ(i)/dVzo dZ(i)/dParam(j) where Z=[X,Y,Z,Vx,Vy,Vz], Param=Force-Parameters | ';
end

% External Orbit comparison
data_product  = 'External Orbit comparison';
test_data_functional = strcmp(data_functional,data_product);
if test_data_functional == 1
    EXT_reference_frame  = 'Orbital Frame';
    test_reference_frame = strcmp(reference_frame,EXT_reference_frame);
    if test_reference_frame == 1    
    format_line = 'Modified Julian Day number | Seconds since 00h | Radial Along-Track Cross-Track | ';
    end
    EXT_reference_frame  = 'ICRF';
    test_reference_frame = strcmp(reference_frame,EXT_reference_frame);
    if test_reference_frame == 1    
    format_line = 'Modified Julian Day number | Seconds since 00h | delta-Z(i) where Z=[X,Y,Z,Vx,Vy,Vz]  | ';
    end
    EXT_reference_frame  = 'Kepler';
    test_reference_frame = strcmp(reference_frame,EXT_reference_frame);
    if test_reference_frame == 1    
    format_line = 'Modified Julian Day number | Seconds since 00h | delta-Z(i) where Z=[a,e,i,Omega,omega,Mean_anomaly] | ';
    end
end

% Intersatellite-Observations
data_product = 'range-rate';
data_functional_keyword = sscanf(data_functional,'%*s %s %*');
test_data_functional = strcmp(data_functional_keyword,data_product);
if test_data_functional == 1
format_line = 'Modified Julian Day number | Seconds since 00h | Range-Rate residuals (m/sec) ';
end

data_product = 'range';
data_functional_keyword = sscanf(data_functional,'%*s %s %*');
test_data_functional = strcmp(data_functional_keyword,data_product);
if test_data_functional == 1
format_line = 'Modified Julian Day number | Seconds since 00h | Range residuals (m)';
end

fprintf(fid,'%-33s %s %s ', 'Data lines format', ': ', format_line);
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% header_end
header_end = 'end_of_header';
fprintf(fid,'%-33s ',header_end );
fprintf(fid,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB orbit format to GNV1b format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO satellite id
test = strcmp(orbiting_object_name,'GRACE-A');
if test == 1
    grace_sat_id = "A";
end
test = strcmp(orbiting_object_name,'GRACE-B');
if test == 1
    grace_sat_id = "B";
end
test = strcmp(orbiting_object_name,'GRACE-C');
if test == 1
    grace_sat_id = "C";
end
test = strcmp(orbiting_object_name,'GRACE-D');
if test == 1
    grace_sat_id = "D";
end

% Reference Frame ID
test = strcmp(reference_frame,'ICRF');
if test == 0
    ref_frame = "E";
elseif test == 1
    ref_frame = "I";
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
data_functional_category = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write computed orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_data_functional = strcmp(data_functional,'Orbit');
if test_data_functional == 1
data_functional_category = 1;
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
fixmjd = 1;
[orbit_matrix_2] = mjd2mjdtt(data_matrix,fixmjd);
[N1, N2] = size(orbit_matrix_2);
Nepochs = N1;
Nelements = N2;
%format_i = ['%' num2str(wno) '.' num2str(prno(i,1)) 'f'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_epochs = 1 : Nepochs
    % Write computed orbit of epoch ti
    data_ith = orbit_matrix_2(i_epochs, :)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
georb_format = 1;
% GEORB orbit data format
if georb_format == 0
    % MJD
    fprintf(fid,'%9d',orbit_matrix_2(i_epochs,1) );
    % Seconds since start of the day  (00h)
    fprintf(fid,'%19.9f',orbit_matrix_2(i_epochs,2) );
    % State Vector: Position Vector
    fprintf(fid,'%29.11f',orbit_matrix_2(i_epochs,3:5) );
    if Nelements > 5
    % State Vector: Velocity Vector
    fprintf(fid,'%29.15f',orbit_matrix_2(i_epochs,6:8) );
    end
    % Chnange line
    fprintf(fid,'%s\n','');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNV1b orbit data format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GEORB orbit format to GNV1b format    
    mjd_TT = (orbit_matrix_2(i_epochs,1));
    mjd_TT = fix(orbit_matrix_2(i_epochs,1));
    sec00_TT = orbit_matrix_2(i_epochs,2);
    % [sec00_UTC,sec00_GPS] = time_scales(sec00_TT,mjd_TT);
    sec00_GPS = sec00_TT - 51.184;
    mjd_gps = mjd_TT;
    % if sec00_GPS < 0
    %     sec00_GPS = 24*3600 + sec00_GPS;
    %     mjd_gps = mjd_gps - 1;
    % end
    % sec00_GPS;
    % mjd_gps;
    % GNV1b orbit data: GPS Time start epoch: 12:00 01-Jan-2000
    [jd2000,mjd2000] = mjd3(12*3600,1,1,2000);
    sec00_GPS_round = round(sec00_GPS);
    time_gnv1b = ( (mjd_gps - mjd2000) * (24*3600) + sec00_GPS_round );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if 1 > 0
    % t_gps2000
    fprintf(fid,'%9d ',time_gnv1b );
    % GRACE-FO sat ID
    fprintf(fid,'%1c ',grace_sat_id );
    % GRACE-FO Reference Frame ID
    fprintf(fid,'%1c ',ref_frame );
    % State Vector: Position Vector
    % fprintf(fid,'%29.11f',orbit_matrix_2(i_epochs,3:5) );    
    %orbit_matrix_2(i_epochs,3)
    fprintf(fid,'%23.11f',orbit_matrix_2(i_epochs,3) );
    %orbit_matrix_2(i_epochs,4)
    fprintf(fid,'%23.11f',orbit_matrix_2(i_epochs,4) );
    %orbit_matrix_2(i_epochs,5)
    fprintf(fid,'%23.11f',orbit_matrix_2(i_epochs,5) );
    % State Vector Errors: Position Vector Errors
    fprintf(fid,'%15.11f', 0,0,0);
    if Nelements > 5
    % State Vector: Velocity Vector
    fprintf(fid,'%25.15f',orbit_matrix_2(i_epochs,6:8) );
    % State Vector Errors: Velocity Vector Errors
    fprintf(fid,'%15.11f', 0,0,0);
    end
    % Flag 10000000
    flag_endline = 10000000;
    fprintf(fid,'% 8d', flag_endline);
    % Chnange line
    fprintf(fid,'%s\n','');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Estimated Parameters matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_data_functional = strcmp(data_functional,'Parameters');
if test_data_functional == 1
data_functional_category = 1;    
    % Data
    fprintf(fid,'%29.15e',data_matrix);
    % Change line
    fprintf(fid,'%s\n','');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data matrix e.g. Partials, Gravity Gradient, Earth Rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_data_functional = strcmp(data_functional,'Partial Derivatives');
% if test_data_functional == 1
% test_data_functional = strcmp(data_functional,'Orbit');
% if test_data_functional == 0  
if data_functional_category == 0
georb_format = 1;
% GEORB orbit data format
if georb_format == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
fixmjd = 1;
[orbit_matrix_2] = mjd2mjdtt(data_matrix,fixmjd);
[N1, N2] = size(orbit_matrix_2);
Nepochs   = N1;
Nelements = N2;
%format_i = ['%' num2str(wno) '.' num2str(prno(i,1)) 'f'];
for i_epochs = 1 : Nepochs
    % Write computed orbit of epoch ti
    data_ith = orbit_matrix_2(i_epochs, :)';
    % MJD
    fprintf(fid,'%9d',orbit_matrix_2(i_epochs,1) );
    % Seconds since start of the day  (00h)
    fprintf(fid,'%19.9f',orbit_matrix_2(i_epochs,2) );    
    % Partial Derivatives
    fprintf(fid,'%29.15e',orbit_matrix_2(i_epochs,3:Nelements) );
    % Chnange line
    fprintf(fid,'%s\n','');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
% GNV1b orbit data format
orbit_matrix_2 = data_matrix;
[N1, N2] = size(orbit_matrix_2);
Nepochs   = N1;
Nelements = N2;
for i_epochs = 1 : Nepochs
    mjd_TT = (orbit_matrix_2(i_epochs,1));
    mjd_TT = fix(orbit_matrix_2(i_epochs,1));
    sec00_TT = orbit_matrix_2(i_epochs,2);
    % [sec00_UTC,sec00_GPS] = time_scales(sec00_TT,mjd_TT);
    sec00_GPS = sec00_TT - 51.184;
    mjd_gps = mjd_TT;
    % if sec00_GPS < 0
    %     sec00_GPS = 24*3600 + sec00_GPS;
    %     mjd_gps = mjd_gps - 1;
    % end
    % sec00_GPS;
    % mjd_gps;
    % GNV1b orbit data: GPS Time start epoch: 12:00 01-Jan-2000
    [jd2000,mjd2000] = mjd3(12*3600,1,1,2000);
    sec00_GPS_round = round(sec00_GPS);
    time_gnv1b = ( (mjd_gps - mjd2000) * (24*3600) + sec00_GPS_round );
    % t_gps2000
    fprintf(fid,'%9d ',time_gnv1b );
    % GRACE-FO sat ID
    fprintf(fid,'%1c ',grace_sat_id );
    % GRACE-FO Reference Frame ID
    fprintf(fid,'%1c ',ref_frame );
    % Partial Derivatives
    fprintf(fid,'%29.15e',orbit_matrix_2(i_epochs,3:Nelements) );   
    % % Flag 10000000
    % flag_endline = 10000000;
    % fprintf(fid,'% 8d', flag_endline);
    % Chnange line
    fprintf(fid,'%s\n',''); 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fclose(fid);
