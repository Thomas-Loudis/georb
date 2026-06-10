function [potential_partials_rpl, potential_partials_xyz] = harmonics_partials_1st(r,n_max,m_max,GM,ae,Cnm,Snm, legendre_functions_struct, n_min) 

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
% Thomas D. Papanikolaou                                     September 2007 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 06/02/2019  Dr. Thomas Papanikolaou 
%             Correction to "dV_r" formula following the loop sum computation   
% 09/12/2022  Thomas Loudis Papanikolaou 
%             Write new function harmonics_partials1.m and minor code upgrade 
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
 
% to = tic; 
 
% 2nd approach: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial derivatives of potential with respect to spherical coordinates : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% - dV_r     : partial derivative of geopotential to radius 
% - dV_phi   : partial derivative of geopotential to latitude 
% - dV_lamda : partial derivative of geopotential to longtitude 
dV_r = 0; 
dV_phi = 0; 
dV_lamda = 0; 
for n = n_min : n_max 
    if n > m_max 
        m_limit = m_max; 
    else 
        m_limit = n; 
    end 
    for m = 0 : m_limit     
        dV_r = dV_r         + ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) ); 
        dV_phi = dV_phi     + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda)); 
        dV_lamda = dV_lamda + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda)); 
    end 
end 
dV_r = (GM/l^2) * dV_r;                 
dV_phi = (GM / l) * dV_phi; 
dV_lamda = (GM / l) * dV_lamda; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
dV_radius = dV_r; 
 
potential_partials_rpl = [dV_radius; dV_phi; dV_lamda]; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial derivatives of (r,phi,lamda) with respect to (x,y,z) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PDVrx = [ 
%       cos(phi)*cos(lamda)            cos(phi)*sin(lamda)          sin(phi) 
%   (-1/l)*sin(phi)*cos(lamda)     (-1/l)*sin(phi)*sin(lamda)    (1/l)*cos(phi) 
% ( -1/(l*cos(phi)) )*sin(lamda)  ( 1/(l*cos(phi)) )*cos(lamda)        0 
% ]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial derivatives of (r,theta,lamda) with respect to (x,y,z) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PDVrx = [ 
%       sin(theta)*cos(lamda)            sin(theta)*sin(lamda)           cos(theta) 
%   ( 1/l)*cos(theta)*cos(lamda)     ( 1/l)*cos(theta)*sin(lamda)    (-1/l)*sin(theta) 
% ( -1/(l*sin(theta)) )*sin(lamda)  ( 1/(l*sin(theta)) )*cos(lamda)         0 
% ]; 
% Replacement of "theta" with "phi" 
PDVrx = [ 
      cos(phi)*cos(lamda)            cos(phi)*sin(lamda)          sin(phi) 
  ( 1/l)*sin(phi)*cos(lamda)     ( 1/l)*sin(phi)*sin(lamda)    (-1/l)*cos(phi) 
( -1/(l*cos(phi)) )*sin(lamda)  ( 1/(l*cos(phi)) )*cos(lamda)        0 
]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Computation of Cartesian counterparts of the acceleration 
fxyz = PDVrx' * [dV_r; dV_phi; dV_lamda]; 
fx = fxyz(1,1); 
fy = fxyz(2,1); 
fz = fxyz(3,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
potential_partials_xyz = fxyz; 
 
% fprintf('%s %.7f \n', 'Time (sec): Original loop:',toc(to)); 
 
% potential_partials_xyz_0 = potential_partials_xyz; 
% potential_partials_rpl_0 = potential_partials_rpl; 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % PAR function degree/order 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% to = tic; 
% n_max_parfor = n_max 
% potential_partials_rpl_degree = zeros(3,n_max+1); 
% potential_partials_xyz_degree = zeros(3,n_max+1); 
% Ncpu_cores = 4 
% parfor (n_par = n_min+1 : n_max+1, Ncpu_cores) 
% degree_par = n_par - 1 
% order_par  = m_max; 
% [potential_partials_rpl_degree(:,n_par), potential_partials_xyz_degree(:,n_par)] = harmonics_partials_1st_degree(r,degree_par,order_par,GM,ae,Cnm,Snm, legendre_functions_struct); 
% end 
% potential_partials_rpl = sum(potential_partials_rpl_degree,2); 
% potential_partials_xyz = sum(potential_partials_xyz_degree,2); 
%  
% fprintf('%s %.7f \n', 'Time (sec): PAR function Ncpu=4:',toc(to)); 
%  
% delta_potential_partials_xyz = potential_partials_xyz - potential_partials_xyz_0 
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
% to = tic; 
% n_max_parfor = n_max 
% potential_partials_rpl_degree = zeros(3,n_max+1); 
% potential_partials_xyz_degree = zeros(3,n_max+1); 
% Ncpu_cores = 0 
% parfor (n_par = n_min+1 : n_max+1, Ncpu_cores) 
% degree_par = n_par - 1 
% order_par  = m_max; 
% [potential_partials_rpl_degree(:,n_par), potential_partials_xyz_degree(:,n_par)] = harmonics_partials_1st_degree(r,degree_par,order_par,GM,ae,Cnm,Snm, legendre_functions_struct); 
% end 
% potential_partials_rpl = sum(potential_partials_rpl_degree,2); 
% potential_partials_xyz = sum(potential_partials_xyz_degree,2); 
%  
% fprintf('%s %.7f \n', 'Time (sec): PAR function Ncpu=0:',toc(to)); 
%  
% delta_potential_partials_xyz = potential_partials_xyz - potential_partials_xyz_0 
%  
%  
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % PAR function order loop 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% to = tic; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Partial derivatives of potential with respect to spherical coordinates : 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % - dV_r     : partial derivative of geopotential to radius 
% % - dV_phi   : partial derivative of geopotential to latitude 
% % - dV_lamda : partial derivative of geopotential to longtitude 
%  
% % dV_r_matrix     = zeros(n_max+1); 
% % dV_phi_matrix   = dV_r_matrix;  
% % dV_lamda_matrix = dV_r_matrix; 
%  
% Ncpu_cores = 4 
% parfor (n_par = n_min+1 : n_max+1, Ncpu_cores) 
%     n = n_par - 1; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     if n > m_max 
%         m_limit = m_max; 
%     else 
%         m_limit = n; 
%     end 
%     % for m = 0 : m_limit     
%     %     dV_r_matrix (n+1,m+1)       = + ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) ); 
%     %     dV_phi_matrix (n+1,m+1)     = + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda)); 
%     %     dV_lamda_matrix (n+1,m+1)   = + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda)); 
%     % end 
%     [dV_r_matrix(1,n_par), dV_phi_matrix(1,n_par), dV_lamda_matrix(1,n_par), potential_partials_rpl_degree(:,n_par)] = harmonics_partials_order(r,n_max,m_limit,GM,ae,Cnm,Snm, legendre_functions_struct)     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% end 
% % [d1,d2] = size(dV_r_matrix) 
% dV_r     = sum(dV_r_matrix,2); 
% dV_phi   = sum(dV_phi_matrix,2); 
% dV_lamda = sum(dV_lamda_matrix,2); 
%  
% dV_r = (GM/l^2) * dV_r;                 
% dV_phi = (GM / l) * dV_phi; 
% dV_lamda = (GM / l) * dV_lamda; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% dV_radius = dV_r; 
%  
% potential_partials_rpl_1 = [dV_radius; dV_phi; dV_lamda]; 
%  
% potential_partials_rpl = sum(potential_partials_rpl_degree,2); 
% potential_partials_rpl = [(GM/l^2) * potential_partials_rpl(1,1); (GM / l) * potential_partials_rpl(2,1); (GM / l) * potential_partials_rpl(3,1)]; 
%  
% fprintf('%s %.7f \n', 'Time (sec): PAR loop Ncpu=4:',toc(to)); 
%  
% delta_potential_partials_rpl = potential_partials_rpl - potential_partials_rpl_0 
%  
% delta_potential_partials_rpl_111 = potential_partials_rpl - potential_partials_rpl_1 
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
% to = tic; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Partial derivatives of potential with respect to spherical coordinates : 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % - dV_r     : partial derivative of geopotential to radius 
% % - dV_phi   : partial derivative of geopotential to latitude 
% % - dV_lamda : partial derivative of geopotential to longtitude 
%  
% dV_r_matrix     = zeros(1,n_max+1); 
% dV_phi_matrix   = dV_r_matrix;  
% dV_lamda_matrix = dV_r_matrix; 
%  
% Ncpu_cores = 0 
% parfor (n_par = n_min+1 : n_max+1, Ncpu_cores) 
%     n = n_par - 1; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     if n > m_max 
%         m_limit = m_max; 
%     else 
%         m_limit = n; 
%     end 
%     [dV_r_matrix(1,n_par), dV_phi_matrix(1,n_par), dV_lamda_matrix(1,n_par), potential_partials_rpl_degree(:,n_par)] = harmonics_partials_order(r,n_max,m_limit,GM,ae,Cnm,Snm, legendre_functions_struct)     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% end 
%  
% dV_r     = sum(dV_r_matrix,2); 
% dV_phi   = sum(dV_phi_matrix,2); 
% dV_lamda = sum(dV_lamda_matrix,2); 
%  
% dV_r = (GM/l^2) * dV_r;                 
% dV_phi = (GM / l) * dV_phi; 
% dV_lamda = (GM / l) * dV_lamda; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% dV_radius = dV_r; 
%  
% fprintf('%s %.7f %s %d \n', 'Time (sec): PAR loop Ncpu:',toc(to),'Ncpu',Ncpu_cores); 
%  
% potential_partials_rpl_1 = [dV_radius; dV_phi; dV_lamda]; 
% potential_partials_rpl = sum(potential_partials_rpl_degree,2); 
% potential_partials_rpl = [(GM/l^2) * potential_partials_rpl(1,1); (GM / l) * potential_partials_rpl(2,1); (GM / l) * potential_partials_rpl(3,1)]; 
%  
%  
% delta_potential_partials_rpl = potential_partials_rpl - potential_partials_rpl_0 
% delta_potential_partials_rpl_111 = potential_partials_rpl - potential_partials_rpl_1 
 
 
