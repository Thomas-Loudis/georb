function [partials_p, partials_c, partials_s] = potential_partials_coef(r,n_max,n_min,GM,ae,gravity_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational potential partials w.r.t. gravitational parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the partial derivatives of the gravitational potential
%  w.r.t. gravitational parameters e.g. spherical harmoncis coefficients
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - r:                  Position vector (m)
% - GM:                 Earth gravity constant  (m^3/sec^2)
% - ae:                 radius  (meters)
% - Cnm, Snm:           normalized spherical harmonics coefficients
% - n_max:              Truncation Degree of harmoncis series expansion 
% - m_max:              Truncation Order of harmoncis series expansion 
%
% Output arguments:
% - partials_c:     partials w.r.t. harmonics coefficients Cnm 
% - partials_s:     partials w.r.t. harmonics coefficients Snm 
% - partials_p:     OVerall matrix of partials w.r.t. harmonics coefficients Cnm and Snm 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Dr. Thomas Loudis Papanikolaou                                     
% Created: 21 March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[partials_p, partials_c, partials_s] = harmonics_partials_coef(r,n_max,n_min,GM,ae,gravity_struct);

