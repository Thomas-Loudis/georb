function [era, gmst] = gmst_era(mjd_TT, mjd_UT)  

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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% gmst_era : GMST and ERA  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - mjd_TT          : MJD in Terrestrial Time including fraction of the day 
% - mjd_UT          : MJD in Universal Time including fraction of the day 
% 
% Output arguments: 
% - era             : Earth Rotation angle in radians 
% - gmst            : Greenwich Mean Sidereal Time (GMST) in radians 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas Loudis Papanikolaou                                  29 July  2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Coefficient for Conversion from arcsec to radians 
radcoef = pi / (180 * 3600); 
const_2PI = 6.283185307179586476925287;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
JD_TT  = mjd_TT + 2400000.5; 
JD_UT1 = mjd_UT + 2400000.5; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Earth Rotation Angle (ERA) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% computaton of ERA in radians 
Tu = JD_UT1 - 2451545.0;
Tu_fraction = mod(Tu, 1);
% era = 2*pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu ); 
era = mod( (const_2PI * (Tu_fraction + 0.7790572732640 + 0.00273781191135448 * Tu )), const_2PI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Greenwich Mean Sidereal Time (GMST) in radians 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% parameter t 
taph = ( JD_TT - 2451545.0 ) / 36525; 
 
% GMST computation in radians 
%gmst_0 = era + 0.014506 * radcoef + 4612.156534 * radcoef * taph + 1.3915817 * radcoef * taph^2 - 0.00000044 * radcoef * taph^3 - 0.000029956 * radcoef * taph^4 - 0.0000000368 * radcoef * taph^5; 
% arcsec
gmst_term1 = 0.014506 + 4612.156534 * taph + 1.3915817 * taph^2 - 0.00000044 * taph^3 - 0.000029956 * taph^4 - 0.0000000368 * taph^5; 
gmst_term1 = mod(gmst_term1, 648000);
gmst_term1_rad = gmst_term1 * radcoef;
% radians
gmst = era + gmst_term1_rad;
% gmst = mod( gmst , 2*pi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
