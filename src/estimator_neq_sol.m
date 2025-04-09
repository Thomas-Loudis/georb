function [Xmatrix, NEQn, NEQu] = estimator_neq_sol(Amatrix, Wmatrix, sigma_obs) 


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
% if n_sigma > 1
% PA_matrix = zeros(n_obs, n_param);
% Pb_matrix = zeros(n_obs, 1);
% for i_obs = 1 : n_obs 
%     w_ij = 1 / (sigma_obs(i_obs,1) )^2;
%     PA_matrix(i_obs,:) = w_ij * Amatrix(i_obs,:);
%     Pb_matrix(i_obs,1) = w_ij * Wmatrix(i_obs,1); 
% end
% 
% % for i_obs = 1 : n_obs 
% %     w_ij = 1 / (sigma_obs(i_obs,1) )^2;
% %     PA_matrix(i_obs,:) = w_ij * Amatrix(i_obs,:);
% %     
% %     for j_param = 1 : n_param
% %         
% %     end
% % end
% end
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
%end
tol2 = 30;
Xmatrix_sol0 = lsqminnorm(NEQn,NEQu,tol2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted solution
if n_sigma > 1
% NEQn = Amatrix' * PA_matrix;
% NEQu = Amatrix' * Pb_matrix;

% [d1_A d2_A] = size(Amatrix)
% [d1_P d2_P] = size(Pmatrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEQn = Amatrix' * Pmatrix * Amatrix;
NEQu = Amatrix' * Pmatrix * Wmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEQn_sol1 = NEQn;
NEQu_sol1 = NEQu;
tol2 = 30;
Xmatrix_sol1 = lsqminnorm(NEQn,NEQu,tol2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N matrix
NEQn = zeros(n_param,n_param);
for n1 = 1 : n_param
    for n2 = 1 : n_param
        Nij_sum = 0;
        for i_obs = 1 : n_obs
           w_ij = 1 / (sigma_obs(i_obs,1) )^2;
           Nij_sum = Amatrix(i_obs,n1) * w_ij * Amatrix(i_obs,n2);
        end
        NEQn(n1,n2) = Nij_sum;
    end
end

% u matrix
NEQu = zeros(n_param,1);
for n1 = 1 : n_param
    u_sum = 0;
    for i_obs = 1 : n_obs
        w_ij = 1 / (sigma_obs(i_obs,1) )^2;
        u_sum = Amatrix(i_obs,n1) * w_ij * Wmatrix(i_obs,1);
    end
    NEQu(n1,1) = u_sum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEQn_sol2 = NEQn;
NEQu_sol2 = NEQu;
tol2 = 30;
Xmatrix_sol2 = lsqminnorm(NEQn,NEQu,tol2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solutions options :: 
NEQn = NEQn_sol1;
NEQu = NEQu_sol1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Equations matrix inverstion algorithm
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix of parameters
% Cx = sigma^2 * inv(NEQn);

% Covariance matrix of errors
% Cv = sigma^2 * (inv(P_matrix) - Amatrix * inv(NEQn) * Amatrix');
% Cv = sigma^2 * ( - Amatrix * inv(NEQn) * Amatrix');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx = 0;
Cv = 0;
