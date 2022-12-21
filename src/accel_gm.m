function [fx,fy,fz] = accel_gm(r,GM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration due to central Earth Gravity field
% Newton Law of Gravity - Keplerian Orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of satellite's acceleration based on Newton's law of gravity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r: position vector in Inertial System (GCRS)
%
% Output arguments:
% - Acceleration in Inertial System (GCRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     September 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 19/12/2022, Thomas Loudis Papanikolaou
%             Gravitational Constant GM has been set as input argument 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% computation of spherical coordinates
[lamda,phi,l] = lamda_phi(r);

% gradient of geopotential V
fr = - GM / (l^2);

% Cartesian counterparts (fx,fy,fz) of acceleration fr
fx = fr * cos(phi)*cos(lamda);
fy = fr * cos(phi)*sin(lamda);
fz = fr * sin(phi);