function [partials_2nd_spher, partials_2nd_xyz, partials_1st_spher, partials_1st_xyz] = potential_partials_2nd(r,n_max,m_max,GM,ae,Cnm,Snm, legendre_functions_struct, n_min)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational potential partials 2nd order w.r.t. position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of the 2nd order partial derivatives of the gravitational
%   potential w.r.t. to spherical and cartestian coordinates
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
% - partials_1st_rpl:   1st order partials w.r.t. spherical coordinates
%                       radius, latitude, longitude
% - partials_1st_xyz:   1st order partials w.r.t. Cartesian coordinates xyz
% - partials_2nd_rpl:   2nd order partials w.r.t. spherical coordinates
%                       radius, latitude, longitude
% - partials_2nd_xyz:   2nd order partials w.r.t. Cartesian coordinates xyz
%                       radius, latitude, longitude
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Dr. Thomas Loudis Papanikolaou                                     
% Created: 9 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[partials_2nd_spher, partials_2nd_xyz, partials_1st_spher, partials_1st_xyz] = harmonics_partials_2nd(r,n_max,m_max,GM,ae,Cnm,Snm, legendre_functions_struct, n_min);

