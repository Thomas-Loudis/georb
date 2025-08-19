function [dV_radius, dV_phi, dV_lamda, potential_partials_rpl] = harmonics_partials_order(r,n_max,m_max,GM,ae,Cnm,Snm, legendre_functions_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational potential partials 1st order w.r.t. position
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
% Thomas Loudis Papanikolaou                                   11 June 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of spherical coordinates in radians
[lamda,phi,l] = lamda_phi(r);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % computation of normalized associated Legendre functions
% [Pnm_norm] = Legendre_functions(phi,n_max);
% % First-order derivatives of normalized associated Legendre functions
% [dPnm_norm] = Legendre1ord(phi,n_max) ;
Pnm_norm = legendre_functions_struct.Pnm_norm;
dPnm_norm = legendre_functions_struct.Pnm_norm_derivatives_1st;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Partial derivatives of potential with respect to spherical coordinates :
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % - dV_r     : partial derivative of geopotential to radius
% % - dV_phi   : partial derivative of geopotential to latitude
% % - dV_lamda : partial derivative of geopotential to longtitude
% dV_r = 0;
% dV_phi = 0;
% dV_lamda = 0;
% for n = n_min : n_max
%     if n > m_max
%         m_limit = m_max;
%     else
%         m_limit = n;
%     end
%     for m = 0 : m_limit    
%         dV_r = dV_r         + ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );
%         dV_phi = dV_phi     + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));
%         dV_lamda = dV_lamda + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));
%     end
% end
% dV_r = (GM/l^2) * dV_r;                
% dV_phi = (GM / l) * dV_phi;
% dV_lamda = (GM / l) * dV_lamda;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dV_radius = dV_r;
% 
% potential_partials_rpl = [dV_radius; dV_phi; dV_lamda];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of potential with respect to spherical coordinates :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - dV_r     : partial derivative of geopotential to radius
% - dV_phi   : partial derivative of geopotential to latitude
% - dV_lamda : partial derivative of geopotential to longtitude
dV_r = 0;
dV_phi = 0;
dV_lamda = 0;

n = n_max;
for m = 0 : m_max
    % dV_r_matrix (n+1,m+1)       = + ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );
    % dV_phi_matrix (n+1,m+1)     = + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));
    % dV_lamda_matrix (n+1,m+1)   = + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));

    dV_r = dV_r         + ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );
    dV_phi = dV_phi     + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));
    dV_lamda = dV_lamda + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));   
end
dV_radius = dV_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
potential_partials_rpl = [dV_radius; dV_phi; dV_lamda];

