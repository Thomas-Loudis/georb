function [orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no, IC_Sec_00, TAI_UTC_table, IAU_PN_XYs_matrix] = prm_ic(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_ic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Configuration file read : 
%  Initial Conditions (State Vector and Initial epoch) and orbit arc length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Orbit configuration file name 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment: This function replaced the function prm1 (01/2011) according to
% the new configuration file format for the orbit modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                                7 July 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2    
% Orbit arc lenght (in hours)
    param_keyword = 'Orbit_arc_length';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    orbit_arc_length_sec = sscanf(param_value,'%f %*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param_keyword = 'Reference_frame';
    [IC_reference_frame] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'Time_scale';
    [IC_time_scale] = read_param_cfg(cfg_fname,param_keyword);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Epoch
    param_keyword = 'Date_format';
    [IC_epoch_date_format] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'Initial_Epoch';
    [param_value, param_line] = read_param_cfg(cfg_fname,param_keyword);
    %IC_epoch_epoch = sscanf(line_ith,'%*s %s %*')
    sec   = sscanf(param_line,'%s %*');
    day   = sscanf(param_line,'%*s %s %*');
    month = sscanf(param_line,'%*s %*s %s %*');
    year  = sscanf(param_line,'%*s %*s %*s %s %*');
    IC_epoch_char = sprintf('%s %s %s %s ',sec,day,month,year);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial State Vector
    param_keyword = 'State_vector';
    [param_value, param_line] = read_param_cfg(cfg_fname,param_keyword);
    IC_Zo(1,1) = sscanf(param_line,'%e %*');
    IC_Zo(2,1) = sscanf(param_line,'%*s %e %*');
    IC_Zo(3,1) = sscanf(param_line,'%*s %*s %e %*');
    IC_Zo(4,1) = sscanf(param_line,'%*s %*s %*s %e %*');
    IC_Zo(5,1) = sscanf(param_line,'%*s %*s %*s %*s %e %*');
    IC_Zo(6,1) = sscanf(param_line,'%*s %*s %*s %*s %*s %e %*');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation modelling | Earth Orientation Parameters (EOP)         
    param_keyword = 'EOP_filename';
    [EOP_filename] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'EOP_interpolation_points';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    EOP_interp_no = sscanf(param_value,'%d %*');

    param_keyword = 'precession_nutation_model';
    [precession_nutation_model] = read_param_cfg(cfg_fname,param_keyword);        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .in configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');
    
% Orbit arc lenght (in hours)
    test = strcmp(str1,'Orbit_arc_length');
    if test == 1
      orbit_arc_length_sec = sscanf(line_ith,'%*s %f %*');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    test = strcmp(str1,'Reference_frame');
    if test == 1
      IC_reference_frame = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Time_scale');
    if test == 1
      IC_time_scale = sscanf(line_ith,'%*s %s %*');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Epoch
    test = strcmp(str1,'Date_format');
    if test == 1
      IC_epoch_date_format = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Initial_Epoch');
    if test == 1
      %IC_epoch_epoch = sscanf(line_ith,'%*s %s %*')
      sec   = sscanf(line_ith,'%*s %s %*');
      day   = sscanf(line_ith,'%*s %*s %s %*');
      month = sscanf(line_ith,'%*s %*s %*s %s %*');
      year  = sscanf(line_ith,'%*s %*s %*s %*s %s %*');
      IC_epoch_char = sprintf('%s %s %s %s ',sec,day,month,year);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial State Vector
    test = strcmp(str1,'State_vector');
    if test == 1
      IC_Zo(1,1) = sscanf(line_ith,'%*s %f %*');
      IC_Zo(2,1) = sscanf(line_ith,'%*s %*s %f %*');
      IC_Zo(3,1) = sscanf(line_ith,'%*s %*s %*s %f %*');
      IC_Zo(4,1) = sscanf(line_ith,'%*s %*s %*s %*s %f %*');
      IC_Zo(5,1) = sscanf(line_ith,'%*s %*s %*s %*s %*s %f %*');
      IC_Zo(6,1) = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %f %*');
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation modelling | Earth Orientation Parameters (EOP)                                                     
    test = strcmp(str1,'EOP_filename');
    if test == 1
      EOP_filename = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'EOP_interpolation_points');
    if test == 1
      EOP_interp_no = sscanf(line_ith,'%*s %d %*');
    end
        
    test = strcmp(str1,'precession_nutation_model');
    if test == 1
      precession_nutation_model = sscanf(line_ith,'%*s %s %*');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit arc length in seconds
orbit_arc_length = orbit_arc_length_sec; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Epoch in Modified Julian Day (MJD) number in TT (Terrestrial Time) time scale 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(IC_epoch_date_format,'MJD');    
if test == 1
    Sec_00 = sscanf(IC_epoch_char,'%*s %*s %f %*');
    MJD_day = sscanf(IC_epoch_char,'%*s %*s %*s %d %*');
    MJDo = MJD_day + Sec_00 / 86400;
end

test = strcmp(IC_epoch_date_format,'calendar');
if test == 1
    IC_sec   = sscanf(IC_epoch_char,'%f %*');
    IC_day   = sscanf(IC_epoch_char,'%*s %d %*');
    IC_month = sscanf(IC_epoch_char,'%*s %*s %d %*');
    IC_year  = sscanf(IC_epoch_char,'%*s %*s %*s %d %*');    
    [JDo, MJDo] = MJD_date(IC_sec, IC_day, IC_month, IC_year);
    Sec_00 = IC_sec;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leap Seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IERS Leap Seconds dat file that provides the Table of TAI-UTC time
% difference for all dates of leap second' introduction by IERS
leap_second_filename = 'Leap_Second.dat';
% Read IERS Leap_Second.dat file
[TAI_UTC, TAI_UTC_table] = read_leapsecond(leap_second_filename, MJDo); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time scale
test = strcmp(IC_time_scale,'GPS');
if test == 1
    dt_TT_GPS_sec = 51.184;
    % dt TT-GPS in days
    dt_TT_GPS = dt_TT_GPS_sec / (24 * 60 * 60);
    % MJDo in TT time scale
    MJDo = MJDo + dt_TT_GPS;
    % Seconds since start of the day (00h) in TT
    Sec_00 = (MJDo - fix(MJDo)) * 86400;     
end

IC_MJDo = MJDo;
IC_Sec_00 = Sec_00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation Parameters data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read EOP data file by IERS C04 Solution
time_period = orbit_arc_length;  % Seconds
[EOP_data, IAU_PN_XYs_matrix] = prm_eop(EOP_filename, EOP_interp_no, time_period, MJDo, precession_nutation_model, TAI_UTC_table);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial State Vector in Intertial reference frame 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(IC_reference_frame,'KEPLER');    
if test == 1
    kepler_0        = IC_Zo;
    M_anomaly_deg   = kepler_0(1,6);
    eccentr         = kepler_0(1,2);
    % Call function kepler_eq for computing Eccentric anomaly
    [E_anomaly_deg] = kepler_eq(M_anomaly_deg,eccentr);
    kepler_0(6)     = E_anomaly_deg;
    % Keplerian elements to state vector Cartesian coordinates
    [zo]            = kepler_z(kepler_0);
    IC_Zo_vec = zo';
else
    IC_Zo_vec = IC_Zo;
end

% Tranformation from ITRF to ICRF
test = strcmp(IC_reference_frame,'ITRF');    
if test == 1
    [eop,deop] = trs2crs(MJDo, EOP_data, EOP_interp_no);
    ro_itrf = [IC_Zo_vec(1,1) IC_Zo_vec(2,1) IC_Zo_vec(3,1)]';
    vo_itrf = [IC_Zo_vec(4,1) IC_Zo_vec(5,1) IC_Zo_vec(6,1)]';
    ro_icrf = eop * ro_itrf;
    vo_icrf = eop * vo_itrf + deop * ro_itrf;
    IC_Zo_vec = [ro_icrf; vo_icrf];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



