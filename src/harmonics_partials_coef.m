function [partials_param, partials_c, partials_s] = harmonics_partials_coef(r,n_max,n_min,GM,ae,gravity_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of the gravitational potential w.r.t. gravity field parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the partial derivatives of the gravitational potential
%  w.r.t. spherical harmonics coefficients
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
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
% Thomas Loudis Papanikolaou                                  22 March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model structure matrix
gfm_struct = gravity_struct;
C_degree_order = gfm_struct.C_degree_order_estim;
S_degree_order = gfm_struct.S_degree_order_estim;
% Gravity Field parameters :: C,S coefficients' degree and order
[Nparam_C, k] = size(C_degree_order);
[Nparam_S, k] = size(S_degree_order);
degree_max = C_degree_order(Nparam_C,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of spherical coordinates in radians
[lamda,phi,l] = lamda_phi(r);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,n_max);
% First-order derivatives of normalized associated Legendre functions
[dPnm_norm] = Legendre1ord(phi,n_max) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% Matrices Pre-allocation
Nparam_Cnm = Nparam_C;
Nparam_Snm = Nparam_S;
partials_c = zeros(3 , Nparam_Cnm);
partials_s = zeros(3 , Nparam_Snm);
partials_param = zeros(3 , Nparam_Cnm+Nparam_Snm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients
% i_coef = 0;
for i_C = 1 : Nparam_C
    n = C_degree_order(i_C, 1);
    m = C_degree_order(i_C, 2);
    
    dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * cos(m*lamda);
    dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ;
    dF_lamda_Cnm =  (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (- sin(m*lamda));    
    % Cartesian counterparts of the partials
    fxyz_Cnm = PDVrx' * [dF_r_Cnm; dF_theta_Cnm; dF_lamda_Cnm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Cnm;
    
    partials_c(:,i_C) = fxyz_Cnm;
end

% S coefficients
for i_C = 1 : Nparam_S
    n = S_degree_order(i_C, 1);
    m = S_degree_order(i_C, 2);

    dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * sin(m*lamda);
    dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
    dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ;    
    % Cartesian counterparts of the partials
    fxyz_Snm = PDVrx' * [dF_r_Snm; dF_theta_Snm; dF_lamda_Snm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Snm;
    partials_s(:,i_C) = fxyz_Snm;
end       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overall matrix for C and S coefficients
partials_param = [partials_c partials_s]; 
