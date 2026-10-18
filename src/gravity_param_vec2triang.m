function [dC_triang,dS_triang] = gravity_param_vec2triang(C_degree_order, S_degree_order, Xmatrix) 

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
% Function: gravity_param_vec2triang 
%-------------------------------------------------------------------------- 
% Purpose: 
%  Gravity Field parameters for harmonics coefficients C,S of degree and
%  order: Converted from vector matrix (Xmatrix) to triangular matrices
%-------------------------------------------------------------------------- 
% Input arguments: 
% - C_degree_order: Degree, Order values of the C coefficients estimated and provided as the initial elements of the given vector matrix Xmatrix
% - S_degree_order: Degree, Order values of the S coefficients estimated and provided as the following elements of the given vector matrix Xmatrix
%  
% Output arguments: 
% - dC_triang:      C coefficients in triangular matrix format
% - dS_triang:      S coefficients in triangular matrix format
% 
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                  18 April 2023 
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Gravity Field parameters :: C,S coefficients' degree and order 
[Nparam_C, k] = size(C_degree_order); 
[Nparam_S, k] = size(S_degree_order); 
degree_max = C_degree_order(Nparam_C,1); 
 
dC_triang = zeros(degree_max+1, degree_max+1);
dS_triang = zeros(degree_max+1, degree_max+1); 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
Nparam = 0; 
 
% Cnm estimated corrections 
for i_C = 1 : Nparam_C 
    Nparam = Nparam + 1; 
    n = C_degree_order(i_C, 1); 
    m = C_degree_order(i_C, 2); 
    dC_triang(n+1,m+1) = Xmatrix(Nparam);
    % dC_triang_i = dC_triang(n+1,m+1)
    % Xmatrix_i = Xmatrix(Nparam)
end 
 
% Snm estimated corrections 
for i_C = 1 : Nparam_S 
    Nparam = Nparam + 1; 
    n = S_degree_order(i_C, 1); 
    m = S_degree_order(i_C, 2); 
    dS_triang(n+1,m+1) = Xmatrix(Nparam); 
end 
%-------------------------------------------------------------------------- 
