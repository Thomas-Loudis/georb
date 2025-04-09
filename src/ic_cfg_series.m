function [ic_data_mjd2] = ic_cfg_series(ic_data_mjd1, i_arc, arc_length)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  ic_cfg_series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read IC configuration file for getting the IC epochs per orbit object id 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - ic_data_mjd1    : Initial Conditions data line (format of IC configuration file)
% - i_arc           : Orbit arc offset multiplier 
% - arc_length      : Orbit arc length in hours 
%
% Output arguments:
% - ic_data_mjd2    : Modified Initial Conditions data line (format of IC configuration file) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             23 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_orbit_arcs = sscanf(ic_data_mjd1,'%*s %*s %*s %*s %*s %*s %*s %*s %d %*');
Orbit_arc_length = arc_length;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test date format
Date_format = sscanf(ic_data_mjd1,'%*s %*s %*s %s %*');

test = strcmp(Date_format,'MJD');
if test == 1
    MJD_day = sscanf(ic_data_mjd1,'%*s%*s%*s%*s %d %*');
    Sec_00h = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s %f %*');
end

test = strcmp(Date_format,'calendar');
if test == 1
    year  = sscanf(ic_data_mjd1,'%*s%*s%*s%*s %d %*');
    month = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s %d %*');
    day   = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s%*s %d %*');
    sec   = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s%*s%*s %f %*');
    [JD_ith, MJD_ith] = MJD_date(sec, day, month, year);
    MJD_day = fix(MJD_ith);
    Sec_00h = sec;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New date MJD
Orbit_arc_length_sec = Orbit_arc_length * 3600;
orbit_arc_days = Orbit_arc_length / 24;

MJD_day = MJD_day + (i_arc-1) * fix(orbit_arc_days);
Sec_00h = Sec_00h + Orbit_arc_length_sec - fix(orbit_arc_days) * 86400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write updated IC data line
object_id  = sscanf(ic_data_mjd1,'%s %*');
ref_frame  = sscanf(ic_data_mjd1,'%*s %s %*');
time_scale = sscanf(ic_data_mjd1,'%*s %*s %s %*');
%IC_Zo = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s%*s%*s%*s%*s %f %f %f %f %f %f %*');
% IC_param = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s%*s%*s%*s%*s %700c %*');
IC_param = sscanf(ic_data_mjd1,'%*s%*s%*s%*s%*s%*s%*s%*s%*s %10000c');

% ic_data_mjd2 = sprintf('%s %s %s %s %d %f 0 0 %d %s ', object_id, ref_frame, time_scale, 'MJD', MJD_day, Sec_00h, N_orbit_arcs, IC_param); 
ic_data_mjd2 = sprintf('%s %s %s %s %d %f %d %d %d %s', object_id, ref_frame, time_scale, 'MJD', MJD_day, Sec_00h, 0, 0, N_orbit_arcs, IC_param) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if i_arc == 1
%     ic_data_mjd2 = ic_data_mjd1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
