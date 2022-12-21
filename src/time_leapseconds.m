function [UTC_TAI] = time_leapseconds(MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Scales Function: time_leapseconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Leap seconds as introduced by IERS bulletin C
%  Provide the time difference between UTC and TAI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - MJD         :   Modified Julian Date of input epoch
%
% Output arguments
% - UTC_TAI     :   Difference between UTC and TAI in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dr. Thomas Papanikolaou                                
% Written: 11 May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/06/2021, Dr. Thomas Papanikolaou
%             Code upgrade
% 07/12/2022, Dr. Thomas Loudis Papanikolaou
%             Code upgrade to read the Leap_Second.dat file provided by
%             IERS as part of the Bulletin C announcements
% 19/12/2022, Dr. Thomas Loudis Papanikolaou
%             Code upgrade: global appraoch to the TAI-UTC table obtained
%             from the Leap_Second.dat file provided by IERS as part of the
%             Bulletin C announcements 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global TAI_UTC_table_glob 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTC (Coordinated Universal Time)
%  Leap seconds are announced by the IERS - Bulletin C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJD to calendar date
[sec,D,M,Y] = MJD_inv(MJD);

TAI_UTC_table = TAI_UTC_table_glob;
[d1 d2] = size(TAI_UTC_table);

for i_step = 1 : d1 
    % MJD number    
    MJD_i = TAI_UTC_table(i_step,1); 
    if MJD >= MJD_i
        % TAI-UTC time difference in seconds    
        TAI_UTC_i = TAI_UTC_table(i_step,2); 
        TAI_UTC = TAI_UTC_i;
    end
end

% UTC-TAI time scales offset
UTC_TAI = - TAI_UTC;
