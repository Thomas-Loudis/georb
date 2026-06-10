function [NEQn, NEQu, Amatrix, bmatrix] = neq_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs) 

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
% Normal Equations for orbit partials  
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
% - orbobs:         Pseudo-Observations based on Kinematic Orbit positions (XYZ)  
%                   orbobs = [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% - sigma_obs       Sigma of observations for forming the weight matrix P 
% 
% Output arguments: 
% - NEQn:           Normal Equations' N matrix  
% - NEQu:           Normal Equations' u matrix  
% - Amatrix:        Design matrix 
% - bmatrix:        b matrix Reduced Observation matrix (Measurements minus Observation function values)   
%  
% Comment:  
% Function neq_orbit is a revised part of the function estimator_orbit  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                             11 August 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Common Epochs of the Observations/Meausurements  
[epochs_common, orbref_common, orbobs_common] = dataseries_commonepochs(orbref,orbobs); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design Matrix :: formulated at the common epochs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
epochs_veq = orbref(:,1); 
[Amatrix] = neq_veq2Amatrix (veqZarray,veqParray, epochs_common, epochs_veq); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% b matrix :: Reduced Observation matrix (Meaurements minus Observation functions values) formed for the common epochs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[bmatrix] = neq_obs2bmatrix (orbref_common,orbobs_common); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Normal Equations matrices 
[NEQn, NEQu] = estimator_neq(Amatrix, bmatrix, sigma_obs); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
