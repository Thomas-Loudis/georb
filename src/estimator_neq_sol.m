function [Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cv] = estimator_neq_sol(Amatrix, Wmatrix, sigma_obs, inv_id) 


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
[NEQn, NEQu] = estimator_neq(Amatrix, Wmatrix, sigma_obs); 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 > 0
% Covariance matrix of parameters
% inv_id = 6; 
[Nmatrix_inv] = inv_mat(NEQn, inv_id);
Cx = sigma^2 * Nmatrix_inv;

% Covariance matrix of errors
% if n_sigma > 1
% Cv  = sigma^2 * (pinv(Pmatrix,Ntol) - Amatrix * pinv(NEQn,Ntol) * Amatrix');
% else
% Cv  = sigma^2 * ( - Amatrix * pinv(NEQn,Ntol) * Amatrix');
% end 
Cv = 0;

else
Cx = 0;
Cv = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
