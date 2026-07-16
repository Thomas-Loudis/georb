function [delta_UT1, delta_LOD] = eop_UT1_libr(mjd_TT, UT1)  

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
% eop_UT1_libr : Libration effect to UT1 and Length of Day (LOD)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - mjd_TT          : MJD in TT including the fraction of the day 
% - UT1             : UT1 time in seconds since the 00h of MJD in TT  
% 
% Output arguments: 
% - delta_UT1       : Libration effect correction to UT1 in seconds 
% - delta_LOD       : Libration effect correction to LOD in sec/day 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas Loudis Papanikolaou                                  29 July  2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Reference: IERS Conventions (2010), Table 5.1b 
% Coefficients for UT1 and LOD libration due to Tidal gravitation   
 
% Collumns   1-6: Arguments Multipliers 
% Collumns     7: Period in days  
% Collumns   8-9: SIN and COS Coefficients for UT1  
% Collumns 10-11: SIN and COS Coefficients for LOD  
%  
% Units: 10^-6 arcsec for Polar motion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tidalgravitation_coefficients_table51b = [ ... 
       2  -2   0  -2   0  -2  0.5377239   0.05  -0.03   -0.3   -0.6  
       2   0   0  -2  -2  -2  0.5363232   0.06  -0.03   -0.4   -0.7  
       2  -1   0  -2   0  -2  0.5274312   0.35  -0.20   -2.4   -4.1  
       2   1   0  -2  -2  -2  0.5260835   0.07  -0.04   -0.5   -0.8  
       2   0   0  -2   0  -1  0.5175645  -0.07   0.04    0.5    0.8  
       2   0   0  -2   0  -2  0.5175251   1.75  -1.01  -12.2  -21.3  
       2   1   0  -2   0  -2  0.5079842  -0.05   0.03    0.3    0.6  
       2   0  -1  -2   2  -2  0.5006854   0.04  -0.03   -0.3   -0.6  
       2   0   0  -2   2  -2  0.5000000   0.76  -0.44   -5.5   -9.6  
       2   0   0   0   0   0  0.4986348   0.21  -0.12   -1.5   -2.6  
       2   0   0   0   0  -1  0.4985982   0.06  -0.04   -0.4   -0.8 ... 
           ]; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
coefficients_table = tidalgravitation_coefficients_table51b; 
[N_tides, n_table_elements] = size(coefficients_table); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
argument_theta_multipliers_table = coefficients_table(:,1); 
arguments_delaunay_multipliers_table = coefficients_table(:,2:6); 
arguments_multipliers_table = coefficients_table(:,1:6); 
 
UT1_sin_terms = coefficients_table(:,8); 
UT1_cos_terms = coefficients_table(:,9); 
 
LOD_sin_terms = coefficients_table(:,10); 
LOD_cos_terms = coefficients_table(:,11); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Delaunay variables (l l' F D Omega) in radians 
[F1, F2, F3, F4, F5] = delaunay_variables(mjd_TT); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Greenwich Mean Sidereal Time (GMST) in radians 
% [era, gmst] = gmst_era(mjd_TT, mjd_UT);
[gmst, era] = gmst_angle(mjd_TT, UT1);

theta_f = mod( (gmst + pi), 2*pi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
delta_UT1 = 0; 
delta_LOD = 0; 
for i_tides = 1 : N_tides 
    % Fundamental Arguments :: Delaunay variables and Multipliers 
    arguments_multipliers = arguments_multipliers_table(i_tides,:); 
    % ksi in radians 
    ksi_tide_ith = arguments_multipliers * [theta_f; F1; F2; F3; F4; F5]; 
    % ksi_tide_ith = arguments_multipliers * [theta_f; -F1; -F2; -F3; -F4; -F5]; 
    ksi_tide_ith = mod(ksi_tide_ith, 2*pi); 
 
    sin_ksi = sin(ksi_tide_ith); 
    cos_ksi = cos(ksi_tide_ith); 
 
    delta_UT1 = delta_UT1 + UT1_sin_terms(i_tides,1) * sin_ksi + UT1_cos_terms(i_tides,1) * cos_ksi; 
    delta_LOD = delta_LOD + LOD_sin_terms(i_tides,1) * sin_ksi + LOD_cos_terms(i_tides,1) * cos_ksi; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Units Conversion 
 
% micro-arcsec to arcsec  
delta_UT1 = delta_UT1 * 10^-6; 
delta_LOD = delta_LOD * 10^-6; 
