function [dC_matrix_aposteriori, dS_matrix_aposteriori, Xaposteriori_P] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, Xmatrix_P)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: gravity_param_cor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Gravity Field parameters update with estimated corrections 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in  
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  18 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters :: C,S coefficients' degree and order
[Nparam_C, k] = size(C_degree_order);
[Nparam_S, k] = size(S_degree_order);
degree_max = C_degree_order(Nparam_C,1);

dC_estim = zeros(degree_max+1, degree_max+1);
dS_estim = zeros(degree_max+1, degree_max+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nparam = 0;
Xaposteriori_P = zeros(Nparam_C + Nparam_S, 1);

% Cnm estimated corrections
for i_C = 1 : Nparam_C
    Nparam = Nparam + 1;
    n = C_degree_order(i_C, 1);
    m = C_degree_order(i_C, 2);
    dC_estim(n+1,m+1) = Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = Xmatrix_P(Nparam) + dC_matrix_apriori(n+1,m+1);                    
end

% Snm estimated corrections
for i_C = 1 : Nparam_S
    Nparam = Nparam + 1;
    n = S_degree_order(i_C, 1);
    m = S_degree_order(i_C, 2);
    dS_estim(n+1,m+1) = Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = Xmatrix_P(Nparam) + dS_matrix_apriori(n+1,m+1);                    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aposteriori matrices of C,S corrections (time variable gravity)
[dC_matrix_aposteriori, dS_matrix_aposteriori] = harmonics_sum(dC_matrix_apriori,dS_matrix_apriori, dC_estim,dS_estim,-1);
