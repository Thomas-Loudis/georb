function [delta_x, delta_y] = eop_pm_libr(mjd_TT, mjd_UT) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eop_pm_libr : Polar Motion Libration effect due to tidal gravitation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mjd             : MJD including the fraction of the day
%
% Output arguments:
% - delta_x         : Libration effect to Polar motion x coordinate in arcsec
% - delta_y         : Libration effect to Polar motion y coordinate in arcsec 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  29 July  2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: IERS Conventions (2010), Table 5.1a
% Coefficients for polar motion libration due to Tidal gravitation  

% Collumns   1-6: Arguments Multipliers
% Collumns     7: Period in days 
% Collumns   8-9: SIN and COS Coefficients for Polar motion x coordinate 
% Collumns 10-11: SIN and COS Coefficients for Polar motion y coordinate 
% 
% Units: 10^-6 arcsec for Polar motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tidalgravitation_coefficients_table51a = [ ...
      0  0 0  0  0 -1 6798.3837    0.0   0.6   -0.1   -0.1
      0 -1 0  1  0  2 6159.1355    1.5   0.0   -0.2    0.1
      0 -1 0  1  0  1 3231.4956  -28.5  -0.2    3.4   -3.9
      0 -1 0  1  0  0 2190.3501   -4.7  -0.1    0.6   -0.9
      0  1 1 -1  0  0 438.35990   -0.7   0.2   -0.2   -0.7
      0  1 1 -1  0 -1 411.80661    1.0   0.3   -0.3    1.0
      0  0 0  1 -1  1 365.24219    1.2   0.2   -0.2    1.4
      0  1 0  1 -2  1 193.55971    1.3   0.4   -0.2    2.9
      0  0 0  1  0  2 27.431826   -0.1  -0.2    0.0   -1.7
      0  0 0  1  0  1 27.321582    0.9   4.0   -0.1   32.4
      0  0 0  1  0  0 27.212221    0.1   0.6    0.0    5.1
      0 -1 0  1  2  1 14.698136    0.0   0.1    0.0    0.6
      0  1 0  1  0  1 13.718786   -0.1   0.3    0.0    2.7
      0  0 0  3  0  3 9.1071941   -0.1   0.1    0.0    0.9
      0  0 0  3  0  2 9.0950103   -0.1   0.1    0.0    0.6
      1 -1 0 -2  0 -1 1.1196992   -0.4   0.3   -0.3   -0.4
      1 -1 0 -2  0 -2 1.1195149   -2.3   1.3   -1.3   -2.3
      1  1 0 -2 -2 -2 1.1134606   -0.4   0.3   -0.3   -0.4
      1  0 0 -2  0 -1 1.0759762   -2.1   1.2   -1.2   -2.1
      1  0 0 -2  0 -2 1.0758059  -11.4   6.5   -6.5  -11.4
      1 -1 0  0  0  0 1.0347187    0.8  -0.5    0.5    0.8
      1  0 0 -2  2 -2 1.0027454   -4.8   2.7   -2.7   -4.8
      1  0 0  0  0  0 0.9972696   14.3  -8.2    8.2   14.3
      1  0 0  0  0 -1 0.9971233    1.9  -1.1    1.1    1.9
      1  1 0  0  0  0 0.9624365    0.8  -0.4    0.4    0.8 ... 
      ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coefficients_table = tidalgravitation_coefficients_table51a;
[N_tides, n_table_elements] = size(coefficients_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
argument_theta_multipliers_table = coefficients_table(:,1);
arguments_delaunay_multipliers_table = coefficients_table(:,2:6);
arguments_multipliers_table = coefficients_table(:,1:6);

xp_sin_terms = coefficients_table(:,8);
xp_cos_terms = coefficients_table(:,9);

yp_sin_terms = coefficients_table(:,10);
yp_cos_terms = coefficients_table(:,11);
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
theta_f = mod( (gmst + pi), 2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_x = 0;
delta_y = 0;
for i_tides = 1 : N_tides
    % Fundamental Arguments :: Delaunay variables and Multipliers
    arguments_multipliers = arguments_multipliers_table(i_tides,:);
    % ksi in radians
    ksi_tide_ith = arguments_multipliers * [theta_f; F1; F2; F3; F4; F5];
    % ksi_tide_ith = arguments_multipliers * [theta_f; -F1; -F2; -F3; -F4; -F5];
    ksi_tide_ith = mod(ksi_tide_ith, 2*pi);

    sin_ksi = sin(ksi_tide_ith);
    cos_ksi = cos(ksi_tide_ith);

    delta_x = delta_x + xp_cos_terms(i_tides,1) * cos_ksi  + xp_sin_terms(i_tides,1) * sin_ksi;
    delta_y = delta_y + yp_cos_terms(i_tides,1) * cos_ksi  + yp_sin_terms(i_tides,1) * sin_ksi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Units Conversion
% micro-arcsec to arcsec 
delta_x = delta_x * 10^-6;
delta_y = delta_y * 10^-6;
