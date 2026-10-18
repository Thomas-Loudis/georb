function [Nzz, Nzg, Ngg, u_z, u_g] = neq_block(NEQ_N, NEQ_u, n_orbparam, n_gravparam)  

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
% Normal Equations' Block matrices   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  NEQ matrices to NEQ block matrices w.r.t. parameters (arc or global
%  parameters)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - NEQn:           Normal Equations N matrix  
% - NEQu:           Normal Equations u matrix  
% 
% Output arguments: 
% - Nzz:            Normal Equations N block matrix w.r.t. orbital elements 
% - Ngg:            Normal Equations N block matrix w.r.t. common parameters e.g. global parameters, gravity parameters 
% - Nzg:            Normal Equations N offdiagonal block matrix (cross-coupling submatrix)
% - u_z:            Normal Equations u block matrix w.r.t. orbital elements
% - u_g:            Normal Equations u block matrix w.r.t. global parameters e.g. gravity parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                              17 April 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NEQ Block matrices :: orbit arc-related parameters and global parameters
% (gravity field parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% N matrix
Nzz = NEQ_N(1:n_orbparam,1:n_orbparam); 
Nzg = NEQ_N(1:n_orbparam,n_orbparam+1:n_orbparam+n_gravparam); 
Ngz = Nzg'; 
% Ngz = NEQ_N(n_orbparam+1:n_orbparam+n_gravparam , 1:n_orbparam);
Ngg = NEQ_N(1+n_orbparam:n_orbparam+n_gravparam,n_orbparam+1:n_orbparam+n_gravparam); 

% u matrix
u_z = NEQ_u(1:n_orbparam,1); 
u_g = NEQ_u(1+n_orbparam:n_orbparam+n_gravparam,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

