function [Emat,Zmat,Kmat] = neq_scaling(Nmat,umat)  

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
% Matrix Scaling / Equilibration for Normal Equations matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Matrix Scaling / Equilibration for reducing condition number to support 
%  Normal Equations matrix inversion   
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - Nmat:          Normal Equations' N matrix  
% - umat:          Normal Equations' u matrix   
% 
% Output arguments: 
% - Emat:          Scaled N matrix 
% - Zmat:          Scaled u matrix   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                               4 April 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%-------------------------------------------------------------------------- 
% Matrix Scaling, Equilibration
%-------------------------------------------------------------------------- 
diag_Nzz = diag(Nmat);
diag_scale_factors = 1.0 ./ sqrt(diag_Nzz);
mat_scaling = diag(diag_scale_factors);
Kmat = mat_scaling;

Emat = Kmat * Nmat * Kmat;
Zmat = Kmat * umat;

%-------------------------------------------------------------------------- 
% asymmetry_error_Emat_pro = max(max(abs(Emat - Emat')))
%-------------------------------------------------------------------------- 
% Symmetry correction 
Emat = (1/2) * (Emat + Emat');
%-------------------------------------------------------------------------- 
% asymmetry_error_Emat_pos = max(max(abs(Emat - Emat')))
