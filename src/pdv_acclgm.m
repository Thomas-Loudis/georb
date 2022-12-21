function [U] = pdv_acclgm(r,GM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of the acceleration due to central Earth Gravity
% field (Keplerian orbit)
%
% Computation of the partial derivatives of the acceleration is based on
% Newton's law of gravity.
%
% This function computes the partial derivatives of the acceleration in the
% Inertial Reference System.
% 
% First the formulas given here yield the partial derivatives of the
% acceleration in an Earth-fixed coordinated system as a function of the
% Earth-fixed position vector.
% Then coordinate transformations are required to obtain the partial
% derivatives in the Inertial coordinate system which is consistent with
% the equation of motion.
%
% The coordinate systems that are used are the Geocentric Celestial
% Refernce System (GCRS) and the Terrestrial Reference System (ITRS).
% Transformation is based on the IERS Conventions 2003 can be described as:
%   [GCRS] = EOP(t) * [ITRS]
%   EOP: Earth Orientation Parameters
%
% Transformation matrix is calculated by function IERS_EOP(TT,D,M,Y).
%
% Acceleration transformation can be calculated as:
%   f_Inertial = EOP(t) * f_Terrestrial
%
% Transformation of the partial derivatives of the acceleration can be
% calculated as:
%   (df/dr)_Inertial = EOP(t) * (df/dr)_Terrestrial * inv( EOP(t) )
%
%
% The partial derivatives of the acceleration f are the second partial
% derivatives of Geopotential U.
% These partial derivatives are noted as:
%   dfx / dy = Uxy,     dfx / dx = Uxx
% They are related 
%   Uxy = Uyx, Uxz = Uzx, Uyz = Uzy and
% Laplace's equation or Laplacian:
%   Uxx + Uyy + Uzz = 0 
%
% Matrix of second partial derivatives is:
%   U = [ Uxx   Uxy   Uxz 
%         Uxy   Uyy   Uyz
%         Uxz   Uyz   Uzz ]  
% 
% Input arguments:
% - t: seconds sinse 0 hours refers to TT (Terrestrial Time)
% - r: position in Inertial System (GCRS)
%
% Output arguments:
% - Partial derivatives of the acceleration in Inertial System (GCRS)
%   are given in the form of second partial derivatives of geopotential
%   Uxx, Uyy, Uzz, Uxy, Uxz, Uyz through the matrix U.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     September 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 29/05/2012  Function's minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computation of spherical coordinates
[lamda,phi,l] = lamda_phi(r);

% coordinates x,y,z
x = r(1,1);
y = r(2,1);
z = r(3,1);

% Computation of second partial derivatives Uxx,Uyy,Uzz,Uxy,Uxz,Uyz
Uxx = 0; Uyy = 0; Uzz = 0; Uxy = 0; Uxz = 0; Uyz = 0;

Uxx = 3 * x^2 - l^2 ;
Uyy = 3 * y^2 - l^2 ;
Uzz = 3 * z^2 - l^2 ;
Uxy = 3 * x * y ;
Uxz = 3 * x * z ;
Uyz = 3 * y * z ;

% Matrix of second partial derivatives of geopotential is:
   Uo = [ Uxx   Uxy   Uxz 
          Uxy   Uyy   Uyz
          Uxz   Uyz   Uzz ] ;   

   U = ( GM / l^5 ) *  Uo ; 
