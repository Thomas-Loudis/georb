function [Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cv] = estimator_neq_sol(Amatrix, Wmatrix, sigma_obs) 


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n_obs n_param] = size(Amatrix);
[n_sigma d2] = size(sigma_obs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weights matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Diagonal matrix 
% w_sigma = 1 / (sigma_obs)^2;
% Pmatrix = w_sigma * eye(n_obs);
% Pmatrix = w_sigma;
if n_sigma > 1
Pmatrix = zeros(n_sigma,n_sigma);
for i_sigma = 1 : n_sigma
    w_ij = 1 / (sigma_obs(i_sigma,1) )^2;
    Pmatrix(i_sigma,i_sigma) = w_ij;
end
%save Pmatrix.out Pmatrix -ASCII -double
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighting matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if n_sigma == 1
% Weighting matrix :: Identiry matrix    
NEQn = Amatrix' * Amatrix;
NEQu = Amatrix' * Wmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted solution
if n_sigma > 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEQn = Amatrix' * Pmatrix * Amatrix;
NEQu = Amatrix' * Pmatrix * Wmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEQn_sol1 = NEQn;
NEQu_sol1 = NEQu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solutions options :: 
NEQn = NEQn_sol1;
NEQu = NEQu_sol1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations matrix inverstion 
% NEQ_matrix_inv = 5;
% Matrix decomposition/factorisation
% 5th approach :: lsqminnorm
tol2 = 30;
Xmatrix5 = lsqminnorm(NEQn,NEQu,tol2);
Xmatrix = Xmatrix5;
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
if 1 < 0
Ntol = rcond(NEQn);

% Covariance matrix of parameters
Cx = sigma^2 * pinv(NEQn, Ntol);

% Covariance matrix of errors
if n_sigma > 1
Cv  = sigma^2 * (pinv(Pmatrix,Ntol) - Amatrix * pinv(NEQn,Ntol) * Amatrix');
else
Cv  = sigma^2 * ( - Amatrix * pinv(NEQn,Ntol) * Amatrix');
end 

else
Cx = 0;
Cv = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
