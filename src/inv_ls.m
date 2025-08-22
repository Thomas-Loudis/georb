function [Xmatrix] = inv_ls(Nmatrix, Umatrix, inv_id) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix Inversion for Least-Squares Solution to Normal Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Matrix Inversion for supporting Least Squares Solution 
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
% Cond_No_inv_ls = cond(Nmatrix);

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
[L,U,P] = lu(NEQn);
y = L \ (P * NEQu);
Xmatrix = U\y;

elseif inv_id == 6
% 6. Cholesky factorisation
% L = chol(NEQn, 'lower');
[L, flag_sym_pos_01] = chol(NEQn, 'lower');
if flag_sym_pos_01 == 0
y = L \ NEQu;
Xmatrix = L' \ y;
else
% LU factorisation
[L,U,P] = lu(NEQn);
y = L \ (P * NEQu);
Xmatrix = U\y;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta_NxU = Nmatrix * Xmatrix - Umatrix;
