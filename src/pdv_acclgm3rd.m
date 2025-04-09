function [pdv3rdbody] = pdv_acclgm3rd(rbody,rsat,GMbody)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of the perturbation due to a 3rd celestial body 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Partial dervatives matrix of the satellites's perturbing acceleration
%  due to celestial bodies (i.e. Sun, Moon, Planets) is computed according
%  to Newton's law of gravity considering the third body as a point mass.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - rbody  : Celestial body's position vector in GCRF (m)
% - rsat   : Satellite's position vector in GCRF (m)
% - GMbody : Celestial body's Gravity constant GM (m^3/sec^2)
%
%   Position vectors refer to the Geocentric Celestial Reference Frame 
%
% Output arguments:
% - pdv3rdbody : Partial derivatives 3x3 matrix for 3rd body (Sun, Moon,
%                Planets) perturbation in GCRF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsat_body : Vector pointing from the Celestial Body to the Satellite 
rsat_body = rsat - rbody;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsat_body vector's distance computation
[lamda_not,phi_not,l_sat_body] = lamda_phi(rsat_body);
clear lamda_not phi_not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Product between vector (rsat-rbody) and its transpose :
% (rsat-rbody) * (rsat-rbody)'
rsat_body_Tproduct = (rsat - rbody) * (rsat - rbody)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives
%I3x3 = eye(3,3);
I3x3 = [
1 0 0
0 1 0
0 0 1];
pdv3rdbody = -GMbody * (1 / l_sat_body^3) * (I3x3 - (3/l_sat_body^2) * rsat_body_Tproduct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
