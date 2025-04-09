function [EOP,dEOP,EOP_inv,dEOP_inv] = trs2crs(mjd,eop,dpint, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: trs2crs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
%  Implementation of the tranformation between Geocentric Celestial
%  Reference System (GCRS) and International Terrestrial Reference System
%  (ITRS) based on the IERS Conventions 2003.
%  The transformation matrices required for the direct and the inverse
%  tranformation are computed for the position and velocity vectors.
%  The derivatives of the transformation matrices are being computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input arguments:
% - mjd:    computation epoch's MJD in Terrestrial Time (TT) scale
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - EOP:       Direct transformation matrix EOP = Q(t)*R(t)*W(t)
% - dEOP:      Derivative of EOP matrix.
%              This matrix is required for the tranformation of velocity
%              vector from ITRS to GCRS.
%              v_CRS = EOP * v_TRS + dEOP * r_TRS
% - EOP_inv:   Inverse transformation matrix (from CRS to TRS)
%              EOP_inv = W(-t) * R(-t) * Q(-t) = EOP ^ T = (EOP)'
% - dEOP_inv:  Derivative of inverse of EOP.
%              This matrix is required for the inverse tranformation of
%              velocity vector from GCRS to ITRS. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Precession, Nutation and Polar motion are considered as constant for the
%  computation of the derivative of EOP and EOP_inv matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
%  Petit, G.; Luzum, B. IERS Conventions 2010, IERS Technical Note No. 36;
%  Verlag des Bundesamts für Kartographie und Geodäsie: Frankfurt am Main,
%  Germany, 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%  10/04/2011  Thomas Papanikolaou
%              Upgrade and function's rename from IERS_EOP.m to iers_rot.m  
%  19/12/2022, Thomas Loudis Papanikolaou
%              Code minor upgrade and rename to trs2crs. The transformation
%              model has considered the IERS Conventions 2010 and its
%              updates
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XYs_IAU200A = orbit_model_struct.IAU_PN_XYs_matrix;
TAI_UTC_table = orbit_model_struct.TAI_UTC_table;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eop_mjd:  Array of values of the EOP days refer in seconds (UTC)
eop_mjd = eop(:,1);
UT1_UTC = eop(:,4);
xpole   = eop(:,2);
ypole   = eop(:,3);
dX_VLBI = eop(:,5);
dY_VLBI = eop(:,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precession-Nutation model IAU2000 (X,Y,s)
X_IAU2000 = XYs_IAU200A(:,2);
Y_IAU2000 = XYs_IAU200A(:,3);
s_CEO = XYs_IAU200A(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOP Interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrangian Interpolation of EOP parameters is performed at UTC time scale
% Civil date (D,Mh,Yr)
[TT,Dy,Mh,Yr] = MJD_inv(mjd);
% computation of UTC time
[UTC,GPS_time] = time_scales(TT,mjd, TAI_UTC_table);
% MJD in UTC time scale
[jd,mjd_int] = MJD_date(UTC,Dy,Mh,Yr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UT1_UTC = UT1 - UTC
[UT1_UTC_int] = interp_Lagrange(eop_mjd,UT1_UTC,mjd_int,dpint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polar motion
[xp_arcsec] = interp_Lagrange(eop_mjd,xpole,mjd_int,dpint);
[yp_arcsec] = interp_Lagrange(eop_mjd,ypole,mjd_int,dpint);
% Conversion from arcsec to radians
xp = (xp_arcsec / 3600) * (pi / 180);
yp = (yp_arcsec / 3600) * (pi / 180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precession-Nutation model IAU2000 or IAU2006/2000A (X,Y,s)
% (X,Y) in radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VLBI corrections (dX,dY) prior to interpolation 
dX = (dX_VLBI / 3600) * (pi / 180);
dY = (dY_VLBI / 3600) * (pi / 180);
% X_IAU2000 = X_IAU2000 + dX;
% Y_IAU2000 = Y_IAU2000 + dY;

% Interpolation 
[Xt_IAU2000] = interp_Lagrange(eop_mjd,X_IAU2000,mjd_int,dpint);
[Yt_IAU2000] = interp_Lagrange(eop_mjd,Y_IAU2000,mjd_int,dpint);
% Xt_IAU2000 = X_IAU2000(2,1);
% Yt_IAU2000 = Y_IAU2000(2,1);

X_PN = Xt_IAU2000;
Y_PN = Yt_IAU2000;

% VLBI corrections (dX,dY) after interpolation 
% [dX_arcsec] = interp_Lagrange(eop_mjd,dX_VLBI,mjd_int,dpint);
% [dY_arcsec] = interp_Lagrange(eop_mjd,dY_VLBI,mjd_int,dpint);
% Interpolation not applied
% dX_arcsec = dX_VLBI(2,1);
% dY_arcsec = dY_VLBI(2,1);
% (dX,dY)EOP in arcsec
% Conversion from arcsec to radians
% dX = (dX_arcsec / 3600) * (pi / 180);
% dY = (dY_arcsec / 3600) * (pi / 180);
% 
% Precession-Nutation model IAU2000 (X,Y,s) + VLBI corrections (dX,dY)
% (X,Y) = (X,Y)IAU2000 + (dX,dY)VLBI
% X_PN = Xt_IAU2000 + dX;
% Y_PN = Yt_IAU2000 + dY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the CEO(Celestial Ephemeris Origin) in the GCRS
% computation of quantity s(t) in radians
[s] = interp_Lagrange(eop_mjd,s_CEO,mjd_int,dpint);
% Interpolation not applied
% s = s_CEO(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transormation matrix computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of UT1 time
% IERS parameter (UT1-UTC) in sec
UT1 = UT1_UTC_int + UTC ;

% computation of Julian Day Number in TT time
[JD_TT,MJD_TT] = MJD_date(TT,Dy,Mh,Yr);

% computation of Julian Day Number in UT1 time
[JD_UT1,MJD_UT1] = MJD_date(UT1,Dy,Mh,Yr);

% parameter t, used in all following expressions (eq.2)
t = ( JD_TT - 2451545.0 ) / 36525;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Rotation Angle
Tu = JD_UT1 - 2451545.0;
%theta = 2* pi * ( 0.7790572732640 + 1.00273781191135448 * Tu );
theta = 2 * pi * ( 0.7790572732640 + Tu + 0.00273781191135448 * Tu );

% Earth Rotation matrix
% computation of  R(t) = R3(-theta)
R3_theta = [ cos(-theta)   sin(-theta)     0
            -sin(-theta)   cos(-theta)     0
                 0             0           1];
             
R_t = R3_theta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the TEO(Terrestrial Ephemeris Origin) in the ITRS
% s' quantity: s' = -47 microarcseconds * t
% converse from microarcseconds to radians
s_TEO = -(47 * 10^(-6) / 3600) * (pi / 180) *  t;  

% computation of matrix R3(-s')
R3_sTEO = [ cos(-s_TEO)  sin(-s_TEO)     0
           -sin(-s_TEO)  cos(-s_TEO)     0
                0             0          1];          
            
% Polar motion xp,yp
% Motion of the CIP (Celestial Intermediate Pole) in the ITRS
% Note: computation of xp,yp including tidal terms is agnored

% computation of matrices R2(xp), R1(yp)
R2_x = [ cos(xp)    0    -sin(xp)
           0        1       0
         sin(xp)    0     cos(xp) ];

R1_y = [   1       0       0
           0    cos(yp)   sin(yp)
           0   -sin(yp)   cos(yp) ];
       
% computation of matrix W(t)
% W(t)=R3(-s')R2(x)R1(y)
W_t = R3_sTEO * R2_x * R1_y; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q(t) matrix
%
% Motion of the CIP (Celestial Intermediate Pole) in the GCRS
% Computation of X, Y
%
% Position of the CEO(Celestial Ephemeris Origin) in the GCRS
% computation of quantity s(t)
% Delaunay variables - fundamental arguments of nutation
% s in radians

% computation of matrix R3(s)
R3_s = [  cos(s)  sin(s)    0
         -sin(s)  cos(s)    0
            0        0      1 ];          

% X,Y,s in radians

% computation of a
a = 1/2 + (1/8)*(X_PN^2+Y_PN^2);

% computation of matrix Q(t)
A_t = [  1-a*X_PN^2      -a*X_PN*Y_PN          X_PN
        -a*X_PN*Y_PN      1-a*Y_PN^2           Y_PN
          -X_PN             -Y_PN       1-a*(X_PN^2+Y_PN^2) ] ;
      
Q_t = A_t * R3_s  ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation matrix from ITRS to GCRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EOP = Q_t * R_t * W_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of derivative of EOP matrix, dEOP = dEOP / dt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = [ 0  -1   0
      1   0   0
      0   0   0 ] ;
% Earth angular velocity w = 0.7292115*10^-4 rad/s Moritz 1980, IAG 1999
omega = 0.7292115 * 10^(-4);
% note: instantaneous computation of dTHETA / dt should replace omega
% THETA = -GAST

dtheta = omega;

% Derivative of EOP matrix: dEOP
%dEOP = dtheta * A_t * P * R3_s * R_t * R3_sTEO * R2_x * R1_y ;
dEOP = dtheta * Q_t * P * R_t * W_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse transformation (from GCRS to ITRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of inverse matrix of EOP (invEOP)
% Inverse matrix of EOP (invEOP)
%invEOP = (W_t ^ -1) * (R_t ^ -1) * (Q_t ^ -1) = (W_t)' * (R_t)' * (Q_t)'
%invEOP = W(-t) * R(-t) * Q(-t) = EOP ^ T = (EOP)'

% Computation of derivative of inverse matrix of EOP (d_invEOP)
% d_invEOP = d_invEOP / dt 
%d_invEOP = omega * W(-t) * P^T * R(-t) * Q(-t) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of inverse of matrix W(t), W(-t) = W_t_inv
% W(-t) = R1(-y)R2(-x)R3(s')

% computation of matrices R2(-xp), R1(-yp)
R2_x_inv = [ cos(-xp)    0    -sin(-xp)
                0       1        0
             sin(-xp)    0     cos(-xp) ] ;

R1_y_inv = [   1       0         0
               0    cos(-yp)   sin(-yp)
               0   -sin(-yp)   cos(-yp) ];

% computation of matrix R3(s')
R3_sTEO_inv = [ cos(s_TEO)  sin(s_TEO)     0
               -sin(s_TEO)  cos(s_TEO)     0
                    0           0          1];          

% computation of matrix W(-t)
W_t_inv = R1_y_inv * R2_x_inv * R3_sTEO_inv ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Rotation matrix
% computation of  R(-t) = R3(theta)
R3_theta_inv = [ cos(theta)   sin(theta)     0
                -sin(theta)   cos(theta)     0
                     0            0          1];           
R_t_inv = R3_theta_inv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of matrix Q(-t)

% computation of matrix R3(-s)
R3_s_inv = [  cos(-s)  sin(-s)    0
             -sin(-s)  cos(-s)    0
                 0        0       1 ];          

% inv(A(t)) = A(-t) = A(-X,-Y)
A_t_inv = [  1-a*X_PN^2     -a*X_PN*Y_PN             -X_PN
            -a*X_PN*Y_PN     1-a*Y_PN^2              -Y_PN
                X_PN            Y_PN          1-a*(X_PN^2+Y_PN^2) ] ;
           
Q_t_inv = R3_s_inv * A_t_inv ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of inverse matrix of EOP (invEOP)
EOP_inv = W_t_inv * R_t_inv * Q_t_inv ;

% Computation of derivative of inverse matrix of EOP (d_invEOP)
%d_invEOP = omega * W(-t) * P^T * R(-t) * Q(-t) ;
Pinv = [  0   1   0
         -1   0   0
          0   0   0 ] ;

dEOP_inv = dtheta * W_t_inv * Pinv * R_t_inv * Q_t_inv ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
