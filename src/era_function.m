function [era] = era_function(mjd_TT, UT1)  

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% era_function : Earth Rotation Angle (ERA)  
%-------------------------------------------------------------------------- 
% Input arguments: 
% - mjd_TT          : MJD in Terrestrial Time including fraction of the day 
% - UT1:            : UT1 in seconds (since mjd_TT 00h) 
% 
% Output arguments: 
% - era             : Earth Rotation angle in radians 
% 
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                  29 July  2025 
%-------------------------------------------------------------------------- 
% Last modified 
% 15/07/2026  Thomas Loudis Papanikolaou 
%             Source Code upgrade for minimizing rotation rounding errors  
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Coefficient for Conversion from arcsec to radians 
radcoef = pi / (180 * 3600); 
const_2PI = 6.283185307179586476925287;
%-------------------------------------------------------------------------- 

% Modified Julian Day Number in UT1 time 
mjd_UT1 = fix(mjd_TT) + UT1 / 86400;

% Julian Day numbers
JD_TT  = mjd_TT  + 2400000.5; 
JD_UT1 = mjd_UT1 + 2400000.5; 
 
%-------------------------------------------------------------------------- 
% Earth Rotation Angle (ERA) 
%-------------------------------------------------------------------------- 
Tu = JD_UT1 - 2451545.0;

Tu_fraction = mod(Tu, 1);
Tu_fraction2 = UT1 / 86400 + 0.5;
Tu_fraction = Tu_fraction2;

% computation of ERA in radians 
% era = 2*pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu ); 
era = mod( (const_2PI * (Tu_fraction + 0.7790572732640 + 0.00273781191135448 * Tu )), const_2PI);
%-------------------------------------------------------------------------- 
