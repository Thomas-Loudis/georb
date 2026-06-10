function [Xmatrix] = inv_ls(Nmatrix, Umatrix, inv_id)  

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
% Matrix Inversion for Least-Squares Solution to Normal Equations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Matrix Inversion for Least-Squares' Normal Equations solution  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - NEQn:           N matrix required for inversion 
% - NEQu:           u matrix where N * X = u  
% - inv_id :        Inversion algorithm ID list:                 
% 1. mldivide or backslash function    
% 2. INV function     
% 3. Moore-Penrose pseudoinverse based on pinv function  
% 4. lsqminnorm algorithm     
% 5. LU factorisation 
% 6. Cholesky factorisation 
% 
% Output arguments: 
% - Xmatrix:        Least Squares method Solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                            15 October 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 18/08/2025  TLP, Code modified and extracted from function estimator_neq_sol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% inv_ls_inv_id = 5 
% inv_id = inv_ls_inv_id; 
% Cond_No_inv_ls = cond(Nmatrix) 
 
NEQn = Nmatrix; 
NEQu = Umatrix; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Least-Squares Solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Normal Equations matrix inversion  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if inv_id == 1 
% 1. mldivide or backslash     
Xmatrix = NEQn \ NEQu;  
 
elseif inv_id == 2 
% 2. INV function     
Xmatrix = inv(NEQn) * NEQu;  
 
elseif inv_id == 3 
% 3. Moore-Penrose pseudoinverse based on pinv function  
Ntol = rcond(NEQn); 
Xmatrix = pinv(NEQn,Ntol) * NEQu; 
 
elseif inv_id == 4 
% 4. lsqminnorm algorithm     
Ntol = rcond(NEQn); 
Xmatrix = lsqminnorm(NEQn,NEQu,Ntol);  
 
elseif inv_id == 5 
% 5. LU factorisation 
[Xmatrix] = inv_lu(Nmatrix, Umatrix); 
 
elseif inv_id == 6 
% 6. Cholesky factorisation 
% L = chol(NEQn, 'lower'); 
[L, flag_sym_pos_01] = chol(NEQn, 'lower'); 
% flag_sym_pos_01 
 
% Symmetric Positive Definite 
if flag_sym_pos_01 == 0 
y = L \ NEQu; 
Xmatrix = L' \ y; 
 
% Non Symmetric Positive Definite 
else 
% LU factorisation 
[Xmatrix] = inv_lu(Nmatrix, Umatrix); 
 
end 
 
elseif inv_id == 7 
% QR factorisation 
[Q,R] = qr(Nmatrix);             
Xmatrix = R \ (Q' * Umatrix);         
 
elseif inv_id == 8 
% Regularisation 
[n1,n2] = size(Nmatrix); 
sigma = 10^-6; 
Nmatrix_regul = Nmatrix + sigma * eye(n1,n2); 
Xmatrix = Nmatrix_regul \ Umatrix; 
 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% delta_NxU = Nmatrix * Xmatrix - Umatrix; 
