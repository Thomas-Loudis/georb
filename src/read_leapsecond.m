function [TAI_UTC, tai_utc_table] = read_leapsecond(leap_second_filename, MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_leapsecond : Read IERS Leap Second dfile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read Leap Seconds data file provided by International Earth Rotation and
%  Reference Systems Service (IERS). 
%  The Leap Seconds data file (Leap_Second.dat) is created by IERS in the
%  frame of the Bulletin C regarding the introduction of leap seconds
% 
% Input arguments
% - leap_second_filename : IERS Leap Seconds data file name that reports
%                          the time difference between TAI and UTC for the
%                          dates of leap second' introduction
% - MJD                  : Modified Julian Date of input epoch
%
% Output arguments:
% - TAI_UTC     : Time difference between TAI and UTC in seconds 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Thomas Loudis Papanikolaou                        
% Created: 7 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Leap Second data file provided by IERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(leap_second_filename);
headeri = 1;
while headeri == 1
    HDlineread = fgetl(fid);
    str_1 = sscanf(HDlineread,'%s%*');
    test = strcmp(str_1,'#');
    if test == 0
        headeri = 0;
    end
end

i_table = 0;
while (~feof(fid))
    lineread = fgetl(fid);
    % MJD number    
    MJD_i = sscanf(lineread,'%f%*');
    % TAI-UTC time difference in seconds    
    TAI_UTC_i = sscanf(lineread,'%*f%*d%*d%*d%d%*');
    if MJD >= MJD_i
        MJD_0 = MJD_i;
        TAI_UTC = TAI_UTC_i;
    end
    % TAI-UTC table values
    i_table = i_table + 1;
    tai_utc_table(i_table,1) = MJD_i;
    tai_utc_table(i_table,2) = TAI_UTC_i;
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
