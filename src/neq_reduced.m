function [NEQ_N_reduced, NEQ_u_reduced, Nzz, Nzg, Ngg, u_z, u_g] = neq_reduced(NEQ_N, NEQ_u, n_orbparam, n_gravparam)  

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
% Reduced Normal Equations' matrices   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Reduced NEQ matrices :: pre-elimination of orbit arc-related parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - NEQn:           Normal Equations N matrix  
% - NEQu:           Normal Equations u matrix  
% 
% Output arguments: 
% - NEQn:           Normal Equations N reduced matrix 
% - NEQu:           Normal Equations u reduced matrix  
% - Xmatrix:        Least Squares method solution  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                              17 April 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 04/09/2026  Thomas Loudis Papanikolaou 
%             Symmetry test and corrections have been added  
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Reduced NEQ matrices :: pre-elimination of orbit arc-related parameters 
%-------------------------------------------------------------------------- 
Nzz = NEQ_N(1:n_orbparam,1:n_orbparam); 
Nzg = NEQ_N(1:n_orbparam,n_orbparam+1:n_orbparam+n_gravparam); 
Ngz = Nzg'; 
% Ngz = NEQ_N(n_orbparam+1:n_orbparam+n_gravparam , 1:n_orbparam);
Ngg = NEQ_N(1+n_orbparam:n_orbparam+n_gravparam,n_orbparam+1:n_orbparam+n_gravparam); 
 
u_z = NEQ_u(1:n_orbparam,1); 
u_g = NEQ_u(1+n_orbparam:n_orbparam+n_gravparam,1); 
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Pro-Scaling
Nzz_0 = Nzz;
Nzg_0 = Nzg;
Ngg_0 = Ngg;
u_z_0 = u_z;
u_g_0 = u_g;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
neq_reduced_NEQscaling = 1;
NEQscaling = neq_reduced_NEQscaling;
if NEQscaling == 1
% Normal Equations scaling
% cond_N_proscale_estimator_neq = cond(Nzz)
% Scaling NEQ matrices N and u
[Nzz,u_z,Kmat] = neq_scaling(Nzz,u_z); 
% cond_N_posscale_estimator_neq = cond(Nzz)
else
   Kmat = 1; 
end

Nzg = Kmat * Nzg;
Ngz = Nzg';
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Symmetry correction 
Nzz = (1/2) * (Nzz + Nzz');

% Cholesky factorisation
inv_mat_id = 6;  
[Nzz_inv] = inv_mat(Nzz, inv_mat_id); 

% Reduction matrix
N_reduction = Ngz * Nzz_inv * Nzg;

% Symmetry correction
N_reduction_sym = (1/2) * (N_reduction + N_reduction');
%-------------------------------------------------------------------------- 

% Pre-elimination of orbital parameters
% NEQ_N_orbelim = (Ngg - Ngz * Nzz_inv * Nzg); 
NEQ_N_orbelim = Ngg - N_reduction_sym; 
NEQ_u_orbelim = u_g - Ngz * Nzz_inv * u_z; 

% Symmetry correction 
NEQ_N_orbelim = (1/2) * (NEQ_N_orbelim + NEQ_N_orbelim');
%-------------------------------------------------------------------------- 

% % Cholesky factorisation
% inv_id = 6; 
% [Xmatrix_grav] = inv_ls(NEQ_N_orbelim, NEQ_u_orbelim, inv_id); 
% Xmatrix = Xmatrix_grav; 
 
NEQ_N_reduced = NEQ_N_orbelim; 
NEQ_u_reduced = NEQ_u_orbelim;

%--------------------------------------------------------------------------
% Output sub-matrices to initial values pro-scaling
Nzz = Nzz_0;
Nzg = Nzg_0;
% Ngg = Ngg_0;
u_z = u_z_0;
% u_g = u_g_0;
%--------------------------------------------------------------------------
