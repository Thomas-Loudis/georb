function [Xmatrix, bmatrix, Amatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cy, Cv] = estimator_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs,inv_id) 

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
% Parameters Estimation 
% Least-squares Estimator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbref:         Dynamic Orbit (observation function)  
%                   orbref = [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% - veqZarray:      Variational Equations solution for the initial state 
%                   vector (State transition matrix) per epoch                    
% - veqParray:      Variational Equations solution for additional parameters 
%                   (force related parameters) per epoch  
% - orbobs:         Pseudo-Observations based on Kinematic Orbit positios (XYZ)  
%                   orbobs = [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% 
% Output arguments: 
% - Xmatrix:        Estimated Parameters corrections  
% - Amatrix:        Design matrix 
% - bmatrix:        b matrix of least squares method   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas D. Papanikolaou                                    August 2013 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 15/04/2021, Thomas Papanikolaou 
%             Upgrade and rename of function mainf_DOD.m 
% 17/08/2022, Thomas Loudis Papanikolaou 
%             Upgrade: matrix decomposition approaches 
% 11/08/2025, TLP,  calling the new function neq_orbit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Formulation of the Normal Equation matrices, Design matrix, b matrix  
[NEQn, NEQu, Amatrix, bmatrix] = neq_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Parameters Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Weighted Least-Squares method solution 
[Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cy, Cv] = estimator_neq_sol(Amatrix, bmatrix, sigma_obs, inv_id); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

