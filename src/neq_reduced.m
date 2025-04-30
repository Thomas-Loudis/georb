function [Xmatrix, NEQ_N_reduced, NEQ_u_reduced] = neq_reduced(NEQ_N, NEQ_u, n_orbparam, n_gravparam) 


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced NEQ matrices :: pre-elimination of orbit arc-related parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nzz = NEQ_N(1:n_orbparam,1:n_orbparam);
Nzg = NEQ_N(1:n_orbparam,n_orbparam+1:n_orbparam+n_gravparam);
Ngz = Nzg';
Ngg = NEQ_N(1+n_orbparam:n_orbparam+n_gravparam,n_orbparam+1:n_orbparam+n_gravparam);

u_z = NEQ_u(1:n_orbparam,1);
u_g = NEQ_u(1+n_orbparam:n_orbparam+n_gravparam,1);

NEQ_N_orbelim = (Ngg - Ngz * pinv(Nzz) * Nzg);
NEQ_u_orbelim =  u_g - Ngz * pinv(Nzz) * u_z;

tol2 = 30;
Xmatrix_grav = lsqminnorm(NEQ_N_orbelim, NEQ_u_orbelim, tol2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NEQ_N_reduced = NEQ_N_orbelim;
NEQ_u_reduced = NEQ_u_orbelim;
Xmatrix = Xmatrix_grav;
