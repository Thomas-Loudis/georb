function [Xmatrix, Xmatrix_alt, bmatrix, Amatrix_out, Cx, Cv, NEQn, NEQu] = estimator_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs,inv_id)


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
[NEQn, NEQu, Amatrix_out, bmatrix] = neq_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cond_No_esitmator_orbit_proscale = cond(NEQn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scaling of Design Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design_matrix_scale_orbit = 0;
if design_matrix_scale_orbit == 1
scale_orbit_pos = 3 * 10^6;
scale_orbit_vel = 3 * 10^3;
[n1,n2] = size(Amatrix_out);

Amatrix_out(:,1:3)  = Amatrix_out(:,1:3) * scale_orbit_pos;
Amatrix_out(:,4:6)  = Amatrix_out(:,4:6) * scale_orbit_vel;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Least-Squares method solution
% [Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cv] = estimator_neq_sol(Amatrix_out, bmatrix, sigma_obs);
[Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cv] = estimator_neq_sol(Amatrix_out, bmatrix, sigma_obs, inv_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cond_No_esitmator_orbit_postscale = cond(NEQn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale parameters
if design_matrix_scale_orbit == 1
Xmatrix(1:3,1) = scale_orbit_pos * Xmatrix(1:3,1);
Xmatrix(4:6,1) = scale_orbit_vel * Xmatrix(4:6,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xmatrix_alt = Xmatrix;

