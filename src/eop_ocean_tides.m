function [delta_x, delta_y, delta_UT1, delta_LOD] = eop_ocean_tides(mjd_TT, mjd_UT)  

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
% eop_ocean_tides : Ocean Tides effect to Polar motion and UT1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - mjd             : MJD including the fraction of the day 
% 
% Output arguments: 
% - delta_x         : Ocean Tidal effect to Polar motion x coordinate in arcsec 
% - delta_y         : Ocean Tidal effect to Polar motion y coordinate in arcsec  
% - delta_UT1       : Ocean Tidal effect correction to UT1 in seconds 
% - delta_LOD       : Ocean Tidal effect correction to LOD in sec/day 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas Loudis Papanikolaou                                  29 July  2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Reference IERS Conventions 2010, Tables 8.2a and 8.2b:  
% Ocean tides' diurnal and semi-dirunal variations coefficients for pole coordinates (xp and yp), UT1 and LOD  
% Collumns   1-6: Arguments Multipliers 
% Collumns   7-8: SIN and COS Coefficients for Polar motion x coordinate  
% Collumns  9-10: SIN and COS Coefficients for Polar motion y coordinate  
% Collumns 11-12: SIN and COS Coefficients for UT1  
%  
% Units: 10^-6 arcse for Polar motion, 10^-6 sec for UT1 
ocean_tides_terms_table8 = [ ... 
      1 -1  0 -2 -2 -2   -0.05    0.94   -0.94   -0.05   0.396  -0.078 ; 
      1 -2  0 -2  0 -1    0.06    0.64   -0.64    0.06   0.195  -0.059 ; 
      1 -2  0 -2  0 -2    0.30    3.42   -3.42    0.30   1.034  -0.314 ; 
      1  0  0 -2 -2 -1    0.08    0.78   -0.78    0.08   0.224  -0.073 ; 
      1  0  0 -2 -2 -2    0.46    4.15   -4.15    0.45   1.187  -0.387 ; 
      1 -1  0 -2  0 -1    1.19    4.96   -4.96    1.19   0.966  -0.474 ; 
      1 -1  0 -2  0 -2    6.24   26.31  -26.31    6.23   5.118  -2.499 ; 
      1  1  0 -2 -2 -1    0.24    0.94   -0.94    0.24   0.172  -0.090 ; 
      1  1  0 -2 -2 -2    1.28    4.99   -4.99    1.28   0.911  -0.475 ; 
      1  0  0 -2  0  0   -0.28   -0.77    0.77   -0.28  -0.093   0.070 ; 
      1  0  0 -2  0 -1    9.22   25.06  -25.06    9.22   3.025  -2.280 ; 
      1  0  0 -2  0 -2   48.82  132.91 -132.90   48.82  16.020 -12.069 ; 
      1 -2  0  0  0  0   -0.32   -0.86    0.86   -0.32  -0.103   0.078 ; 
      1  0  0  0 -2  0   -0.66   -1.72    1.72   -0.66  -0.194   0.154 ; 
      1 -1  0 -2  2 -2   -0.42   -0.92    0.92   -0.42  -0.083   0.074 ; 
      1  1  0 -2  0 -1   -0.30   -0.64    0.64   -0.30  -0.057   0.050 ; 
      1  1  0 -2  0 -2   -1.61   -3.46    3.46   -1.61  -0.308   0.271 ; 
      1 -1  0  0  0  0   -4.48   -9.61    9.61   -4.48  -0.856   0.751 ; 
      1 -1  0  0  0 -1   -0.90   -1.93    1.93   -0.90  -0.172   0.151 ; 
      1  1  0  0 -2  0   -0.86   -1.81    1.81   -0.86  -0.161   0.137 ; 
      1  0 -1 -2  2 -2    1.54    3.03   -3.03    1.54   0.315  -0.189 ; 
      1  0  0 -2  2 -1   -0.29   -0.58    0.58   -0.29  -0.062   0.035 ; 
      1  0  0 -2  2 -2   26.13   51.25  -51.25   26.13   5.512  -3.095 ; 
      1  0  1 -2  2 -2   -0.22   -0.42    0.42   -0.22  -0.047   0.025 ; 
      1  0 -1  0  0  0   -0.61   -1.20    1.20   -0.61  -0.134   0.070 ; 
      1  0  0  0  0  1    1.54    3.00   -3.00    1.54   0.348  -0.171 ; 
      1  0  0  0  0  0  -77.48 -151.74  151.74  -77.48 -17.620   8.548 ; 
      1  0  0  0  0 -1  -10.52  -20.56   20.56  -10.52  -2.392   1.159 ; 
      1  0  0  0  0 -2    0.23    0.44   -0.44    0.23   0.052  -0.025 ; 
      1  0  1  0  0  0   -0.61   -1.19    1.19   -0.61  -0.144   0.065 ; 
      1  0  0  2 -2  2   -1.09   -2.11    2.11   -1.09  -0.267   0.111 ; 
      1 -1  0  0  2  0   -0.69   -1.43    1.43   -0.69  -0.288   0.043 ; 
      1  1  0  0  0  0   -3.46   -7.28    7.28   -3.46  -1.610   0.187 ; 
      1  1  0  0  0 -1   -0.69   -1.44    1.44   -0.69  -0.320   0.037 ; 
      1  0  0  0  2  0   -0.37   -1.06    1.06   -0.37  -0.407  -0.005 ; 
      1  2  0  0  0  0   -0.17   -0.51    0.51   -0.17  -0.213  -0.005 ; 
      1  0  0  2  0  2   -1.10   -3.42    3.42   -1.09  -1.436  -0.037 ; 
      1  0  0  2  0  1   -0.70   -2.19    2.19   -0.70  -0.921  -0.023 ; 
      1  0  0  2  0  0   -0.15   -0.46    0.46   -0.15  -0.193  -0.005 ; 
      1  1  0  2  0  2   -0.03   -0.59    0.59   -0.03  -0.396  -0.024 ; 
      1  1  0  2  0  1   -0.02   -0.38    0.38   -0.02  -0.253  -0.015 ; 
      2 -3  0 -2  0 -2   -0.49   -0.04    0.63    0.24  -0.089  -0.011 ; 
      2 -1  0 -2 -2 -2   -1.33   -0.17    1.53    0.68  -0.224  -0.032 ; 
      2 -2  0 -2  0 -2   -6.08   -1.61    3.13    3.35  -0.637  -0.177 ; 
      2  0  0 -2 -2 -2   -7.59   -2.05    3.44    4.23  -0.745  -0.222 ; 
      2  0  1 -2 -2 -2   -0.52   -0.14    0.22    0.29  -0.049  -0.015 ; 
      2 -1 -1 -2  0 -2    0.47    0.11   -0.10   -0.27   0.033   0.013 ; 
      2 -1  0 -2  0 -1    2.12    0.49   -0.41   -1.23   0.141   0.058 ; 
      2 -1  0 -2  0 -2  -56.87  -12.93   11.15   32.88  -3.795  -1.556 ; 
      2 -1  1 -2  0 -2   -0.54   -0.12    0.10    0.31  -0.035  -0.015 ; 
      2  1  0 -2 -2 -2  -11.01   -2.40    1.89    6.41  -0.698  -0.298 ; 
      2  1  1 -2 -2 -2   -0.51   -0.11    0.08    0.30  -0.032  -0.014 ; 
      2 -2  0 -2  2 -2    0.98    0.11   -0.11   -0.58   0.050   0.022 ; 
      2  0 -1 -2  0 -2    1.13    0.11   -0.13   -0.67   0.056   0.025 ; 
      2  0  0 -2  0 -1   12.32    1.00   -1.41   -7.31   0.605   0.266 ; 
      2  0  0 -2  0 -2 -330.15  -26.96   37.58  195.92 -16.195  -7.140 ; 
      2  0  1 -2  0 -2   -1.01   -0.07    0.11    0.60  -0.049  -0.021 ; 
      2 -1  0 -2  2 -2    2.47   -0.28   -0.44   -1.48   0.111   0.034 ; 
      2  1  0 -2  0 -2    9.40   -1.44   -1.88   -5.65   0.425   0.117 ; 
      2 -1  0  0  0  0   -2.35    0.37    0.47    1.41  -0.106  -0.029 ; 
      2 -1  0  0  0 -1   -1.04    0.17    0.21    0.62  -0.047  -0.013 ; 
      2  0 -1 -2  2 -2   -8.51    3.50    3.29    5.11  -0.437  -0.019 ; 
      2  0  0 -2  2 -2 -144.13   63.56   59.23   86.56  -7.547  -0.159 ; 
      2  0  1 -2  2 -2    1.19   -0.56   -0.52   -0.72   0.064   0.000 ; 
      2  0  0  0  0  1    0.49   -0.25   -0.23   -0.29   0.027  -0.001 ; 
      2  0  0  0  0  0  -38.48   19.14   17.72   23.11  -2.104   0.041 ; 
      2  0  0  0  0 -1  -11.44    5.75    5.32    6.87  -0.627   0.015 ; 
      2  0  0  0  0 -2   -1.24    0.63    0.58    0.75  -0.068   0.002 ; 
      2  1  0  0  0  0   -1.77    1.79    1.71    1.04  -0.146   0.037 ; 
      2  1  0  0  0 -1   -0.77    0.78    0.75    0.45  -0.064   0.017 ; 
      2  0  0  2  0  2   -0.33    0.62    0.65    0.19  -0.049   0.018 ] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
[N_tides, n_table_elements] = size(ocean_tides_terms_table8); 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
argument_theta_multipliers_table = ocean_tides_terms_table8(:,1); 
arguments_delaunay_multipliers_table = ocean_tides_terms_table8(:,2:6); 
arguments_multipliers_table = ocean_tides_terms_table8(:,1:6); 
 
xp_sin_terms = ocean_tides_terms_table8(:,7); 
xp_cos_terms = ocean_tides_terms_table8(:,8); 
 
yp_sin_terms = ocean_tides_terms_table8(:,9); 
yp_cos_terms = ocean_tides_terms_table8(:,10); 
 
UT1_sin_terms = ocean_tides_terms_table8(:,11); 
UT1_cos_terms = ocean_tides_terms_table8(:,12); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Coefficient for Conversion from arcsec to radians 
radcoef = pi / (180 * 3600); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Delaunay variables (l l' F D Omega) in radians 
[F1, F2, F3, F4, F5] = delaunay_variables(mjd_TT); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Greenwich Mean Sidereal Time (GMST) in radians 
% [thetag] = iers_gmst(mjd,eop,dpint, orbit_model_struct); 
[era, gmst] = gmst_era(mjd_TT, mjd_UT); 
 
% Theta angle 
theta_f = mod( (gmst + pi), 2*pi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
delta_x = 0; 
delta_y = 0; 
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
 
% Polar Motion corrections due to Ocean Tides variations 
    delta_x = delta_x + xp_cos_terms(i_tides,1) * cos_ksi  + xp_sin_terms(i_tides,1) * sin_ksi; 
    delta_y = delta_y + yp_cos_terms(i_tides,1) * cos_ksi  + yp_sin_terms(i_tides,1) * sin_ksi; 
 
    % EOP Corrections due to ocean tides terms 
    delta_UT1 = delta_UT1 + UT1_cos_terms(i_tides,1) * cos_ksi  + UT1_sin_terms(i_tides,1) * sin_ksi; 
end 
 
% Units Conversion 
 
% arcsec  
delta_x = delta_x * 10^-6; 
delta_y = delta_y * 10^-6; 
 
% 10^-6 s to seconds  
delta_UT1 = delta_UT1 * 10^-6; 
 
% % 10^-5 s to seconds  
% delta_LOD = delta_LOD * 10^-5; 
