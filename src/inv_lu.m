function [Xmatrix] = inv_lu(Nmatrix, Umatrix) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix Inversion for Least-Squares Solution to Normal Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Matrix Inversion based on LU factorisation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - NEQn:           N matrix required for inversion
% - NEQu:           u matrix where N * X = u 
%
% Output arguments:
% - Xmatrix:        Least Squares method Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             26 August 2025
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
% 5. LU factorisation
LU_scale = 0;

if LU_scale == 0
% LU factorisation with permutation
[L,U,P] = lu(NEQn);
y = L \ (P * NEQu);
Xmatrix = U\y;

elseif LU_scale == 1 
% LU factorisation for sparse matrices; Permutation and diagonal Scaling
[L,U,P,Q,D] = lu(NEQn);
% Scaling 
u_scaled = D * NEQu;
% Permuation 
P_perm = P * u_scaled;
y = L \ P_perm;
Xmatrix = Q * (U \ y);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta_NxU = Nmatrix * Xmatrix - Umatrix;
