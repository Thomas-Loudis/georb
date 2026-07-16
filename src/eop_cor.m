function [delta_UT1, delta_xp, delta_yp] = eop_cor(mjd_TT, UT1) 

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
% Function: eop_cor 
%-------------------------------------------------------------------------- 
% Purpose 
%  EOP corrections due to Tides effects (zonal tides, ocean tides,
%  libration, tidal gravitation)  
%-------------------------------------------------------------------------- 
% 
% Input arguments: 
% - mjd_TT:        MJD in Terrestrial Time (TT) scale 
% - UT1:           UT1 time in seconds since the 00h of MJD in TT (mjd_TT) 
% 
% Output arguments: 
% - delta_UT1:     Corrections to UT1 in seconds
% - delta_xp:      Corrections to polar motion x coordinate in radians
% - delta_yp:      Corrections to polar motion y coordinate in radians
% 
%--------------------------------------------------------------------------  
% Thomas Loudis Papanikolaou                                   10 July 2026 
%-------------------------------------------------------------------------- 
 

% Coefficient for Conversion from arcsec to radians 
const_2PI = 6.283185307179586476925287;
arcsec2rad = const_2PI /1296000; 

%-------------------------------------------------------------------------- 
% EOP variations due to Zonal Tides terms :: applied only for special studies  
% delta_UT1_zonaltides is applied as remove/restore approach during interpolation of the EOP Data :: Moved to eop_interp.m function    
[delta_UT1_zonaltides, delta_LOD_zonaltides, delta_omega_zonaltides] = eop_zonal_tides(mjd_TT); 
%-------------------------------------------------------------------------- 
delta_omega = delta_omega_zonaltides; 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% EOP corrections due to Ocean Tides terms (arcsec, sec) 
[delta_x_oceantides, delta_y_oceantides, delta_UT1_oceantides, delta_LOD_oceantides] = eop_ocean_tides(mjd_TT, UT1); 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% UT1 and LOD libration effect in sec, sec/day for UT1, LOD respectively 
[delta_UT1_libr, delta_LOD_libr] = eop_UT1_libr(mjd_TT, UT1); 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Polar Motion libration correction due to tidal gravitation (arcsec) 
[delta_x_libr, delta_y_libr] = eop_pm_libr(mjd_TT, UT1);  
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% UT1 corrections at interpolation epoch 
delta_UT1   = delta_UT1_oceantides + delta_UT1_libr; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Polar Motion corrections
%-------------------------------------------------------------------------- 
delta_xp_arcsec = delta_x_oceantides + delta_x_libr; 
delta_yp_arcsec = delta_y_oceantides + delta_y_libr; 
 
% Conversion from arcsec to radians 
delta_xp_rad = (delta_xp_arcsec * arcsec2rad); 
delta_yp_rad = (delta_yp_arcsec * arcsec2rad); 

delta_xp = delta_xp_rad;
delta_yp = delta_yp_rad;
%-------------------------------------------------------------------------- 

