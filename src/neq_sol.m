function [Xmatrix, NEQn, NEQu, NEQscaling_mat] = neq_sol(NEQn, NEQu, inv_id, condtol)  

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
% - NEQn:           Normal Equations N matrix  
% - NEQu:           Normal Equations u matrix  
% - inv_id:         Matrix inversion algorithm ID number  
% - condtol:        Matrix condition number tolerance for applying matrix scaling prior inversion  
% 
% Output arguments: 
% - Xmatrix:        Least Squares method solution (estimated parameters correction to apriori values)  
% - NEQn:           Normal Equations N matrix (scaled in case NEQ scaling applied) 
% - NEQu:           Normal Equations u matrix (scaled in case NEQ scaling applied) 
% - NEQscaling_mat: Normal Equations scaling matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                          18 September 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 

%-------------------------------------------------------------------------- 
% Normal Equations scaling
%-------------------------------------------------------------------------- 
ncond_Nmatrix = cond(NEQn);
% ncond_N_proscale_function_neq_sol = cond(NEQn)
NEQscaling_01 = 0;
NEQscaling_ncond = condtol;

if ncond_Nmatrix > NEQscaling_ncond
% Normal Equations scaling
[NEQn, NEQu, Kmat] = neq_scaling(NEQn,NEQu);
% cond_N_postscale_Function_neq_sol = cond(NEQn)
NEQscaling_01 = 1;
else
Kmat = 1;    
end
NEQscaling_mat = Kmat;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Normal Equations Solution for Parameter Estimation 
%-------------------------------------------------------------------------- 
% Normal Equations matrix inversion  
% inv_id = 6; 
[Xmatrix] = inv_ls(NEQn, NEQu, inv_id); 
%-------------------------------------------------------------------------- 

if NEQscaling_01 == 1
% de-scaling of estimated Parameter corrections 
Xmatrix = Kmat * Xmatrix; 
end
