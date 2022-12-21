function [gmst] = iers_gmst(mjd,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: iers_gmst.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
%  Computation of the Greenwich Mean Sidereal Time (GMST) based on the 
%  IAU 2006/2000A precesiion-nutaion model according to the IERS Conventio-
%  ns 2010.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mjd   : computation epoch's MJD in Terrestrial Time (TT) scale
% - eop   : Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint : Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - GMST  : Greenwich Mean Sidereal Time (GMST) in radians
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                    June 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOP : UT1_UTC Interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eop_mjd:  Array of values of the EOP days refer in seconds (UTC)
eop_mjd = eop(:,1);
UT1_UTC = eop(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrangian Interpolation of EOP parameters is performed at UTC time scale
% Civil date (D,Mh,Yr)
[TT,Dy,Mh,Yr] = MJD_inv(mjd);
% computation of UTC time
[UTC,GPS_time] = time_scales(TT,mjd);
% MJD in UTC time scale
[jd,mjd_int] = MJD_date(UTC,Dy,Mh,Yr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UT1_UTC = UT1 - UTC
[UT1_UTC_int] = interp_Lagrange(eop_mjd,UT1_UTC,mjd_int,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of UT1 time
% IERS parameter (UT1-UTC) in sec
UT1 = UT1_UTC_int + UTC ;
% computation of Julian Day Number in UT1 time
[JD_UT1,MJD_UT1] = MJD_date(UT1,Dy,Mh,Yr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Rotation Angle (ERA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computaton of ERA in radians
Tu = JD_UT1 - 2451545.0;
era = 2 * pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient for Conversion from arcsec to radians
radcoef = pi / (180 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greenwich Mean Sidereal Time (GMST) computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of Julian Day Number in TT time
[JD_TT,MJD_TT] = MJD_date(TT,Dy,Mh,Yr);
% parameter t
taph = ( JD_TT - 2451545.0 ) / 36525;
% GMST computation in radians
gmst = era + 0.014506 * radcoef + 4612.156534 * radcoef * taph + 1.3915817 * radcoef * taph^2 - 0.00000044 * radcoef * taph^3 - 0.000029956 * radcoef * taph^4 - 0.0000000368 * radcoef * taph^5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
