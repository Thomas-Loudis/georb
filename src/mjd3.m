function [JD,MJD] = mjd3(t,D,M,Y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion between the Modified Julian Day Number (MJD) and the Calendar
% Date.
%
% Conversion from Civil date to Julian day number and MJD.
% Civil Date is expressed in Gregorian calendar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - t: time in seconds since 0 hours
% - D: number of day
% - M: number of month
% - Y: Year
%
% Output arguments:
% - JD:  Julian day number
% - MJD: Modified Julian day number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converse from given epoch t to UT
% UT: time in hours, fraction of the day
UT = t / 3600;

% computation of Julian day number (JD)
JD = 367*Y - fix( (7/4)*(Y+ fix( (M+9)/12)) ) - fix( (3/4)*(fix((Y+(M-9)/7)/100)+1) ) + fix(275 * M/9) + D + 1721028.5 + (UT/24) ;

% computation of Modified Julian Day Number (MJD)
MJD = JD - 2400000.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
