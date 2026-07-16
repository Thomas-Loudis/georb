function [trs2crs_mat,dtrs2crs_mat, crs2trs_mat,dcrs2trs_mat] = trs2crs_matrix(mjd_TT, UT1, X_PN, Y_PN, s_PN, xp, yp) 

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
% Function: trs2crs_matrices 
%-------------------------------------------------------------------------- 
% Purpose 
%  Tranformation between Geocentric Celestial Reference System (GCRS) and
%  International Terrestrial Reference System (ITRS) based on the IERS
%  Conventions 2010.
%  The transformation matrices required for the direct and the inverse 
%  tranformation are computed for the position and velocity vectors. 
%  The derivatives of the transformation matrices are being computed. 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - mjd_TT:       computation epoch's MJD in Terrestrial Time (TT) scale 
% - UT1:          computation epoch's seconds (since MJD 00h) in UT1 time scale 
% - X_PN:         Precession-Nutation model X coordinate 
% - Y_PN:         Precession-Nutation model Y coordinate 
% - s_PN:         Precession-Nutation model CIO' s coordinate 
% - xp, yp:       Polar motion coordinates 
% 
% Output arguments: 
% - trs2crs_mat:  ITRS to GCRS transformation matrix Q(t)*R(t)*W(t) 
% - dtrs2crs_mat: Derivative of ITRS to GCRS transformation matrix
% - crs2trs_mat:  GCRS to ITRS transformation matrix  W(-t) * R(-t) * Q(-t) 
% - dcrs2trs_mat: Derivative of GCRS to ITRS transformation matrix
%-------------------------------------------------------------------------- 
% Remark: 
%  Precession, Nutation and Polar motion are considered as constant for the 
%  computation of the derivative of EOP and EOP_inv matrices. 
%-------------------------------------------------------------------------- 
% Thomas D. Papanikolaou                                      November 2007 
%-------------------------------------------------------------------------- 
% Last modified: 
%  10/07/2026  Thomas Loudis Papanikolaou 
%              Code upgraded from trs2crs function 
%-------------------------------------------------------------------------- 
 

% Coefficient for Conversion from arcsec to radians 
const_2PI = 6.283185307179586476925287;
arcsec2rad = const_2PI /1296000; 

s = s_PN;

%-------------------------------------------------------------------------- 
% t parameter 
%-------------------------------------------------------------------------- 
% computation of Julian Day Number in TT time 
% [JD_TT,MJD_TT] = MJD_date(TT,Dy,Mh,Yr); 
JD_TT = mjd_TT + 2400000.5; 

% basic parameter t used in Equations (eq.2) 
taph = ( JD_TT - 2451545.0 ) / 36525; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Earth Rotation Angle 
%-------------------------------------------------------------------------- 
[era] = era_function(mjd_TT, UT1);
theta = era;
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Earth Rotation matrix 
%-------------------------------------------------------------------------- 
% computation of  R(t) = R3(-theta)
[R3_theta] = r3_mat(-theta);
R_t = R3_theta; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Position of the TEO (Terrestrial Ephemeris Origin) in the ITRS 
%-------------------------------------------------------------------------- 
% s' quantity: s' = -47 microarcseconds * t 
% converse from microarcseconds to radians 
s_TEO = -47 * 10^(-6) *  taph * arcsec2rad;   
% s_TEO = -1 * mod( (47 * 10^(-6) / 3600) * (pi / 180) * taph , 2*pi);   
 
% computation of matrix R3(-s') 
[R3_sTEO] = r3_mat(-s_TEO);
             
% Polar motion xp,yp 
% Motion of the CIP (Celestial Intermediate Pole) in the ITRS  
% computation of matrices R2(xp), R1(yp) 
[R2_x] = r2_mat(xp);
[R1_y] = r1_mat(yp);
        
% computation of matrix W(t) 
% W(t)=R3(-s')R2(x)R1(y) 
W_t = R3_sTEO * R2_x * R1_y;  
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Q(t) matrix 
%--------------------------------------------------------------------------% 
% Motion of the CIP (Celestial Intermediate Pole) in the GCRS 
% 
% Position of the CEO(Celestial Ephemeris Origin) in the GCRS  
% computation of matrix R3(s) 
[R3_s] = r3_mat(s);
 
% X,Y,s in radians 
 
% computation of a 
a = 1/2 + (1/8)*(X_PN^2+Y_PN^2); 
 
% computation of matrix Q(t) 
A_t = [  1-a*X_PN^2      -a*X_PN*Y_PN          X_PN 
        -a*X_PN*Y_PN      1-a*Y_PN^2           Y_PN 
          -X_PN             -Y_PN       1-a*(X_PN^2+Y_PN^2) ] ; 
       
Q_t = A_t * R3_s ;  
%-------------------------------------------------------------------------- 
if 1 > 0
R2 = X_PN * X_PN + Y_PN * Y_PN ;
if (R2 > 0) 
    E_angle = arctan(Y_PN, X_PN) ;
else
    E_angle = 0;
end
d_angle = atan( sqrt ( R2 / (1-R2) ) ) ;

% Q(t) = R3(-E) * R2(-d) * R3(E) * R3(s);

% E angle R3 matrix
[R3_mE] = r3_mat(-E_angle);
[R3_E] = r3_mat(E_angle);

% D angle R2 matrix
[R2_mD] = r2_mat(-d_angle);
[R2_D] = r2_mat(d_angle);

Q_t2 = R3_mE * R2_mD * R3_E * R3_s;
% delta_Q_t = Q_t2 - Q_t;
Q_t = Q_t2;
end
%-------------------------------------------------------------------------- 
 

%-------------------------------------------------------------------------- 
% Transformation matrix from ITRS to GCRS 
%-------------------------------------------------------------------------- 
EOP = Q_t * R_t * W_t; 
%-------------------------------------------------------------------------- 
 
%-------------------------------------------------------------------------- 
% Computation of derivative of EOP matrix, dEOP = dEOP / dt  
%-------------------------------------------------------------------------- 
P = [ 0  -1   0 
      1   0   0 
      0   0   0 ] ;

% Earth angular velocity w = 0.7292115*10^-4 rad/s Moritz 1980, IAG 1999 
% omega = 0.7292115 * 10^(-4); 
% dtheta = omega; 
dtheta = 2*pi * 1.00273781191135448 * (1/86400) ; 
 
% Zonal Tides Correction 
% dtheta = dtheta + delta_omega; 
 
% Derivative of EOP matrix: dEOP 
% dEOP = dtheta * A_t * P * R3_s * R_t * R3_sTEO * R2_x * R1_y ; 
dEOP = dtheta * Q_t * P * R_t * W_t; 
%-------------------------------------------------------------------------- 
 
 
%-------------------------------------------------------------------------- 
% Inverse transformation (from GCRS to ITRS) 
%-------------------------------------------------------------------------- 
% Computation of inverse matrix of EOP (invEOP) 
% Inverse matrix of EOP (invEOP) 
% invEOP = (W_t^-1) * (R_t^-1) * (Q_t^-1) = (W_t)' * (R_t)' * (Q_t)' 
% invEOP = W(-t) * R(-t) * Q(-t) = EOP^T = (EOP)' 
 
% Computation of derivative of inverse matrix of EOP (d_invEOP) 
% d_invEOP = invEOP / dt  
% d_invEOP = omega * W(-t) * P^T * R(-t) * Q(-t) ; 
%-------------------------------------------------------------------------- 

% Computation of inverse of matrix W(t), W(-t) = W_t_inv 
% W(-t) = R1(-y)R2(-x)R3(s') 
 
% computation of matrices R2(-xp), R1(-yp) 
[R2_x_inv] = r2_mat(-xp);
[R1_y_inv] = r1_mat(-yp);
 
% computation of matrix R3(s') 
[R3_sTEO_inv] = r3_mat(s_TEO);
 
% computation of matrix W(-t) 
W_t_inv = R1_y_inv * R2_x_inv * R3_sTEO_inv ; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Earth Rotation matrix 
% computation of  R(-t) = R3(theta) 
[R3_theta_inv] = r3_mat(theta);
R_t_inv = R3_theta_inv; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% computation of matrix Q(-t) 
 
% computation of matrix R3(-s) 
[R3_s_inv] = r3_mat(-s);
 
% inv(A(t)) = A(-t) = A(-X,-Y) 
A_t_inv = [  1-a*X_PN^2     -a*X_PN*Y_PN             -X_PN 
            -a*X_PN*Y_PN     1-a*Y_PN^2              -Y_PN 
                X_PN            Y_PN          1-a*(X_PN^2+Y_PN^2) ] ; 
            
Q_t_inv = R3_s_inv * A_t_inv ;

% Q_t_inv = R3_s_inv * R3_mE * R2_D * R3_E ;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Computation of inverse matrix of EOP (invEOP) 
EOP_inv = W_t_inv * R_t_inv * Q_t_inv ; 
 
% Computation of derivative of inverse matrix of EOP (d_invEOP) 
%d_invEOP = omega * W(-t) * P^T * R(-t) * Q(-t) ; 
Pinv = [  0   1   0 
         -1   0   0 
          0   0   0 ] ; 
 
dEOP_inv = dtheta * W_t_inv * Pinv * R_t_inv * Q_t_inv ; 
%-------------------------------------------------------------------------- 

% Final output arguments assignment
trs2crs_mat  = EOP;
dtrs2crs_mat = dEOP;
crs2trs_mat  = EOP_inv;
dcrs2trs_mat = dEOP_inv;
