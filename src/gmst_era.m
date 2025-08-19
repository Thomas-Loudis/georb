function [era, gmst] = gmst_era(mjd_TT, mjd_UT) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmst_era : GMST and ERA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mjd             : MJD including the fraction of the day
%
% Output arguments:
% - delta_UT1       : Libration effect correction to UT1 in seconds
% - delta_LOD       : Libration effect correction to LOD in sec/day
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  29 July  2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient for Conversion from arcsec to radians
radcoef = pi / (180 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JD_TT  = mjd_TT + 2400000.5;
JD_UT1 = mjd_UT + 2400000.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Rotation Angle (ERA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computaton of ERA in radians
Tu = JD_UT1 - 2451545.0;
% era = 2 * pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu );
era = mod( (2*pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu )), 2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greenwich Mean Sidereal Time (GMST) in radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [thetag] = iers_gmst(mjd,eop,dpint, orbit_model_struct);

% parameter t
taph = ( JD_TT - 2451545.0 ) / 36525;

% GMST computation in radians
gmst_0 = era + 0.014506 * radcoef + 4612.156534 * radcoef * taph + 1.3915817 * radcoef * taph^2 - 0.00000044 * radcoef * taph^3 - 0.000029956 * radcoef * taph^4 - 0.0000000368 * radcoef * taph^5;
gmst = mod(gmst_0 , 2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IERS Code Formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TWOPI   = 6.283185307179586476925287D0;
RAD2SEC = 86400D0/TWOPI;
sec2rad = 1 / RAD2SEC;

RMJD  = mjd_TT;
RMJD0 = 51544.5D0;                
T = (RMJD-RMJD0)/36525D0;

      GMST = mod (   67310.54841D0 + ...
                    T*( (8640184.812866D0 + 3155760000D0) + ...
                    T*( 0.093104D0 + ...
                    T*( -0.0000062 ))), 86400D0 );


      gmst_iers = GMST / RAD2SEC; % + PI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gmst = gmst_iers;

% delta_GMST_sec = gmst - gmst_iers