function [MJD_clock, calendar_clock] = time_filename()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: time_filename.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Time argument based on computer clock to be used for output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                    6 July 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Current epoch (computer time)
[clock_time] = clock;

% Calendar date
year  = clock_time(1);
month = clock_time(2);
day   = clock_time(3);
sec   = clock_time(4) * 3600 + clock_time(5) * 60 + clock_time(6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified Julian Day number MJD
% Calendar to MJD day number
[JD_clock,MJD_clock] = mjd3(sec,day,month,year);
%MJD_clock_char = sprintf('%d%s%d',fix(MJD_clock),'.',mjd_fr); 
MJD_clock_char = sprintf('%.5f',MJD_clock); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time name based on clock calendar date
clock_char = sprintf('%d',clock_time(1));
for i = 2 : 5
if clock_time(i) < 10 
    clock_i = sprintf('%s%d','0',clock_time(i));
else
    clock_i = sprintf('%d',clock_time(i));
end
clock_char = sprintf('%s%s',clock_char,clock_i);
end    

sec = fix(clock_time(6));
if sec < 10 
    clock_i = sprintf('%s%d','0',sec);
else
    clock_i = sprintf('%d',sec);
end
clock_char = sprintf('%s%s',clock_char,clock_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final output argument for end time
MJD_clock = MJD_clock_char;
calendar_clock = clock_char;
