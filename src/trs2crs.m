function [trs2crs_mat,dtrs2crs_mat, crs2trs_mat,dcrs2trs_mat] = trs2crs(mjd,eop,dpint, orbit_model_struct) 

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
% Function: trs2crs 
%-------------------------------------------------------------------------- 
% Purpose 
%  Implementation of the tranformation between Geocentric Celestial 
%  Reference System (GCRS) and International Terrestrial Reference System 
%  (ITRS) based on the IERS Conventions 2003. 
%  The transformation matrices required for the direct and the inverse 
%  tranformation are computed for the position and velocity vectors. 
%  The derivatives of the transformation matrices are being computed. 
% 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - mjd:    computation epoch's MJD in Terrestrial Time (TT) scale 
% - eop:    Earth Orientation Parameters (EOP) data that are required for 
%           the orbit arc length 
% - dpint:  Number of data points (days) that are required for the EOP 
%           interpolation to the computation epoch 
% 
% Output arguments:
% - trs2crs_mat:  ITRS to GCRS transformation matrix Q(t)*R(t)*W(t) 
% - dtrs2crs_mat: Derivative of ITRS to GCRS transformation matrix
% - crs2trs_mat:  GCRS to ITRS transformation matrix  W(-t) * R(-t) * Q(-t) 
% - dcrs2trs_mat: Derivative of GCRS to ITRS transformation matrix
% 
%-------------------------------------------------------------------------- 
% References: 
%  Petit, G.; Luzum, B. IERS Conventions 2010, IERS Technical Note No. 36; 
%  Verlag des Bundesamts für Kartographie und Geodäsie: Frankfurt am Main, 
%  Germany, 2010. 
%--------------------------------------------------------------------------  
% Thomas D. Papanikolaou, AUTH                                November 2007 
%--------------------------------------------------------------------------  
% Last modified: 
%  10/04/2011  Thomas Papanikolaou 
%              Upgrade and function's rename from IERS_EOP.m to iers_rot.m   
%  19/12/2022, Thomas Loudis Papanikolaou 
%              Code minor upgrade and rename to trs2crs. The transformation 
%              model has considered the IERS Conventions 2010 and its 
%              updates 
% 07/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
% 30/06/2026  Thomas Loudis Papanikolaou 
%             Precesssion-Nutation model minor changes via new code   
% 10/07/2026  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
%-------------------------------------------------------------------------- 
 
 
TAI_UTC_table = orbit_model_struct.TAI_UTC_table; 

%-------------------------------------------------------------------------- 
% Time scales
%-------------------------------------------------------------------------- 
% Civil date (D,Mh,Yr) 
[TT,Dy,Mh,Yr] = MJD_inv(mjd); 

% computation of UTC time 
[UTC,GPS_time] = time_scales(TT,mjd, TAI_UTC_table);

% MJD in UTC time scale 
% [jd,mjd_int] = MJD_date(UTC,Dy,Mh,Yr); 
mjd_int = fix(mjd) + UTC / 86400;
mjd_UTC = mjd_int; 
mjd_TT = mjd; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% EOP Interpolation based on Lagrangian Interpolator
%-------------------------------------------------------------------------- 
% Interpolation of EOP parameters in UTC time scale 
[xp, yp, UT1_UTC_int, dX, dY] = eop_interp(eop, dpint, mjd_UTC) ;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% UT1 time 
%-------------------------------------------------------------------------- 
% (UT1-UTC) in sec 
UT1 = UT1_UTC_int + UTC;
% computation of Julian Day Number in UT1 time 
% [JD_UT1,MJD_UT1] = MJD_date(UT1,Dy,Mh,Yr);
MJD_UT1 = fix(mjd) + UT1 / 86400;
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Precession-Nutation model IAU2006/2000A (X,Y,s) 
%-------------------------------------------------------------------------- 
% X,Y,s in radians by IAU2006/2000A Precession-Nutation model
% Position of the CIO (Celestial Intermediate Origin) in the GCRS 
[X_pn_model, Y_pn_model, s_pn_model] = pn_iau2006_2000a(mjd_TT);
s_PN = s_pn_model;

% XYs_IAU200A = orbit_model_struct.IAU_PN_XYs_matrix; 
% PN_model_function = 1;
% if PN_model_function == 0
% % Precession-Nutation model tables 
% X_IAU2000 = XYs_IAU200A(:,2); 
% Y_IAU2000 = XYs_IAU200A(:,3); 
% s_CEO = XYs_IAU200A(:,4); 
% % Interpolation  
% [X_pn_model] = interp_Lagrange(eop_mjd,X_IAU2000,mjd_int,dpint); 
% [Y_pn_model] = interp_Lagrange(eop_mjd,Y_IAU2000,mjd_int,dpint); 
% [s_PN] = interp_Lagrange(eop_mjd,s_CEO,mjd_int,dpint);
% end
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Precession-Nutation model (X,Y,s) + VLBI corrections (dX,dY) 
X_PN = X_pn_model + dX; 
Y_PN = Y_pn_model + dY; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Tidal Variations to EOP  
%-------------------------------------------------------------------------- 
delta_UT1 = 0; 
delta_xp = 0; 
delta_yp = 0;  
EOP_cor_01 = 1; 
if EOP_cor_01 > 0     
% mjd_UT = mjd_UTC;
mjd_UT = MJD_UT1;
% EOP corrections :: UT1 in sec, Polar Motion (xp,yp) in radians
[delta_UT1, delta_xp, delta_yp] = eop_cor(mjd_TT, mjd_UT);
end 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Polar Motion coordinates (corrected at interpolation epoch) 
xp = xp + delta_xp; 
yp = yp + delta_yp; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% UT1 :: Tidal Variations corrections to UT1 time scale
UT1 = UT1 + delta_UT1; 
% Modified Julian Day Number in UT1 time 
% [JD_UT1,MJD_UT1] = MJD_date(UT1,Dy,Mh,Yr); 
MJD_UT1 = fix(mjd) + UT1 / 86400;
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Transformation matrix computations 
%--------------------------------------------------------------------------
% Transformation matrix from ITRS to GCRS :: trs2crs_mat matrix and dtrs2crs_mat is its derivative w.r.t. time 
% Transformation matrix from GCRS to ITRS :: crs2trs_mat matrix and its derivative dcrs2trs_mat
[trs2crs_mat,dtrs2crs_mat,crs2trs_mat,dcrs2trs_mat] = trs2crs_matrix(mjd_TT, UT1, X_PN, Y_PN, s_PN, xp, yp);
%-------------------------------------------------------------------------- 
