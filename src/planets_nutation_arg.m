function [F6,F7,F8,F9,F10,F11,F12,F13,F14] = planets_nutation_arg(MJD) 

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
% Function: planets_nutation_arg 
%-------------------------------------------------------------------------- 
% Purpose 
%  Computation of the arguments for the planetary nutation based on the
%  IERS Conventions 2010 
% 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - mjd             : MJD including the fraction of the day 
% 
% Output arguments: 
% - F6 to F14        : Planetary Nutation Arguments in radians 
% 
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                   30 June 2026 
%-------------------------------------------------------------------------- 

 
%-------------------------------------------------------------------------- 
% Coefficient for Conversion from arcsec to radians 
arcsec2rad = pi / (180 * 3600); 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% parameter t in Julian Centuries
JD_TT = MJD + 2400000.5; 
taph = ( JD_TT - 2451545.0 ) / 36525; 
%-------------------------------------------------------------------------- 
 

%-------------------------------------------------------------------------- 
% Planetary Nutation Arguments in radians; t in Julian Centuries
%-------------------------------------------------------------------------- 
% F6   LMe 
F6 = 4.402608842 + 2608.7903141574 * taph;

% F7   LV e 
F7 = 3.176146697 + 1021.3285546211 * taph;

% F8   LE 
F8 = 1.753470314 + 628.3075849991 * taph;

% F9   LMa 
F9 = 6.203480913 + 334.0612426700 * taph;

% F10   LJ 
F10 = 0.599546497 + 52.9690962641 * taph;

% F11   LSa 
F11 = 0.874016757 + 21.3299104960 * taph;

% F12   LU 
F12 = 5.481293872 + 7.4781598567 * taph;

% F13   LNe 
F13 = 5.311886287 + 3.8133035638 * taph;

% F14   pA 
F14 = 0.02438175 * taph + 0.00000538691 * taph^2;
%-------------------------------------------------------------------------- 

