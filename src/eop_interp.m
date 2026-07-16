function [xp, yp, UT1_UTC_int, dX, dY] = eop_interp(eop_data, dpint, mjd) 

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
% Function: eop_interp 
%-------------------------------------------------------------------------- 
% Purpose 
%  EOP data interpolation 
%-------------------------------------------------------------------------- 
% 
% Input arguments: 
% - mjd:           MJD of interpolation epoch 
% - eop_data:      EOP data array; format 
%                  eop_data = [ MJD_UTC xpole ypole UT1_UTC dX_VLBI dY_VLBI ]  
% 
% Output arguments: 
% - xp:            polar motion x coordinate in radians
% - yp:            polar motion y coordinate in radians
% - UT1_UTC_int:   UT1-UTC time difference 
% - dX, dY:        VLBI Corrections to the Precession-Nutation model 
% 
%--------------------------------------------------------------------------  
% Thomas Loudis Papanikolaou                                   10 July 2026 
%-------------------------------------------------------------------------- 
 

% Coefficient for Conversion from arcsec to radians 
const_2PI = 6.283185307179586476925287;
arcsec2rad = const_2PI /1296000; 

%-------------------------------------------------------------------------- 
% eop_mjd:  Array of values of the EOP days refer in (UTC) 
%-------------------------------------------------------------------------- 
eop = eop_data;

eop_mjd = eop(:,1); 
UT1_UTC = eop(:,4); 
xpole   = eop(:,2); 
ypole   = eop(:,3); 
dX_VLBI = eop(:,5); 
dY_VLBI = eop(:,6); 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% EOP Interpolation 
%-------------------------------------------------------------------------- 
% Lagrangian Interpolation of EOP parameters is performed at UTC time scale 
mjd_int = mjd;
 
%-------------------------------------------------------------------------- 
% UT1_UTC = UT1 - UTC 
%-------------------------------------------------------------------------- 
EOP_zonaltides = 1;
if EOP_zonaltides == 0
    [UT1_UTC_int] = interp_Lagrange(eop_mjd,UT1_UTC,mjd_int,dpint); 

elseif EOP_zonaltides == 1 
    [n_mjd, n2] = size(eop_mjd); 
    [n_UT1_UTC, n4] = size(UT1_UTC); 
    for i_mjd = 1 : n_mjd 
        mjd_day_data = eop_mjd(i_mjd,1); 
        % EOP variations due to Zonal Tides terms   
        % delta_UT1_zonaltides is removed and then restored to the EOP Data  
        [delta_UT1_zonaltides, delta_LOD_zonaltides, delta_omega_zonaltides] = eop_zonal_tides(mjd_day_data); 
        % delta_UT1   = - delta_UT1_zonaltides + delta_UT1_oceantides; 
        % Remove Zonal TIdes 
        UT1_UTC_ZonalTides_removed(i_mjd,1) = UT1_UTC(i_mjd,1) - delta_UT1_zonaltides; 
    end  
    % Interpolate at compuation epoch     
    [UT1_UTC_int_ZonalTides_no] = interp_Lagrange(eop_mjd,UT1_UTC_ZonalTides_removed,mjd_int,dpint); 
    
    % Restore Zonal Tides effect 
    [delta_UT1_zonaltides, delta_LOD_zonaltides, delta_omega_zonaltides] = eop_zonal_tides(mjd_int);  
    UT1_UTC_int = UT1_UTC_int_ZonalTides_no + delta_UT1_zonaltides; 
end 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Polar motion 
[xp_arcsec] = interp_Lagrange(eop_mjd,xpole,mjd_int,dpint); 
[yp_arcsec] = interp_Lagrange(eop_mjd,ypole,mjd_int,dpint); 
% Conversion from arcsec to radians 
xp = xp_arcsec * arcsec2rad ; 
yp = yp_arcsec * arcsec2rad ;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% VLBI corrections (dX,dY) after interpolation  
[dX_arcsec] = interp_Lagrange(eop_mjd,dX_VLBI,mjd_int,dpint); 
[dY_arcsec] = interp_Lagrange(eop_mjd,dY_VLBI,mjd_int,dpint); 
% Conversion from arcsec to radians 
dX = dX_arcsec * arcsec2rad ; 
dY = dY_arcsec * arcsec2rad ; 
%-------------------------------------------------------------------------- 
