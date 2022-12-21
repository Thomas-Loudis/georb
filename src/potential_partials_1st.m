function [partials_rpl, partials_xyz] = potential_partials_1st(r,n_max,m_max,GM,ae,Cnm,Snm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational potential partials 1st order w.r.t. position vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of the 1st order partial derivatives of the gravitational
%   potential w.r.t. to spherical and cartestian coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r:                  Position vector (m) r = [x y z]'
% - GM:                 Earth gravity constant  (m^3/sec^2)
% - ae:                 radius  (meters)
% - Cnm, Snm:           normalized spherical harmonics coefficients
% - n_max:              Truncation Degree of harmoncis series expansion 
% - m_max:              Truncation Order of harmoncis series expansion 
%
% Output arguments:
% - potential_partials_rpl: 1st order partials w.r.t. spherical coordinates 
%                           radius, latitude, longitude
%   potential_partials_rpl = [partial_radius; partial_latitude; partial_longitude]
%
% - potential_partials_xyz: 1st order partials w.r.t. Cartesian coordinates 
%   potential_partials_xyz = [partial_x; partial_y; partial_z]
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Dr. Thomas Loudis Papanikolaou                                     
% Created: 9 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[partials_rpl, partials_xyz] = harmonics_partials_1st(r,n_max,m_max,GM,ae,Cnm,Snm);

