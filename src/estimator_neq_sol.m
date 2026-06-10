function [Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cy, Cv] = estimator_neq_sol(Amatrix, Wmatrix, sigma_obs, inv_id)  

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
% Least-Squares solution to Normal Equations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Normal Equations and its solution based on Least-Squares method  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - Amatrix:        Design matrix (partials matrix) 
% - Wmatrix:        b matrix (reduced observations matrix)   
% - sigma_obs:      Observations variances matrix (collumn matrix); Variances per observation obtained from the full covariance matrix 
% 
% Output arguments: 
% - Xmatrix:        Least Squares method solution (estimated parameters correction to apriori values)  
% - NEQn:           Normal Equations N matrix  
% - NEQu:           Normal Equations u matrix  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                            15 October 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 15/10/2022  Dr. Thomas Loudis Papanikolaou 
%             Upgrade according to the function estimator_orbit_intersat.m   
% 17/02/2023  Dr. Thomas Loudis Papanikolaou 
%             Upgrade to the weigthts of the NEQ solution   
% 18/08/2025  TLP, Code modification 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Normal Equations 
[NEQn, NEQu, Pmatrix] = estimator_neq(Amatrix, Wmatrix, sigma_obs);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Parameter Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Normal Equations matrix inversion  
% inv_id = 6; 
[Xmatrix] = inv_ls(NEQn, NEQu, inv_id); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
bmatrix = Wmatrix; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Errors matrix 
error_matrix = bmatrix - Amatrix * Xmatrix; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Covariance matrix of estimated parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Number of Observations 
[dim1 dim2] = size(bmatrix); 
n_obs = dim1; 
 
% Number of Parameters  
[dim1 dim2] = size(Xmatrix); 
m_param = dim1; 
 
% Reference Variance 
sigma = sqrt(error_matrix' * error_matrix / (n_obs - m_param) ); 
sigma0 = sigma; 
 
% Covariance matrix of parameters 
% inv_id = 6;  
[Nmatrix_inv] = inv_mat(NEQn, inv_id); 
Cx_matrix = sigma^2 * Nmatrix_inv; 
 
% Covariance matrix of observations (measurements) 
[n_sigma_obs,n2] = size(sigma_obs); 
if n_sigma_obs == 1 
Cy_matrix = sigma^2 * (sigma_obs(1,1) )^2; 
else 
Cy_matrix = sigma^2 * inv_mat(Pmatrix, inv_id); 
end 
 
% Covariance matrix of Observation functions 
% Cyf_matrix = sigma^2 * Amatrix * Nmatrix_inv * Amatrix'; 
 
% Covariance matrix of errors 
% Cv_matrix  = sigma^2 * ( inv_mat(Pmatrix, inv_id) - Amatrix * Nmatrix_inv * Amatrix'); 
Cv_matrix = 0; 
% end  
 
% Diagonal Elements only 
Cx = diag(Cx_matrix); 
Cy = diag(Cy_matrix); 
Cv = diag(Cv_matrix); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
