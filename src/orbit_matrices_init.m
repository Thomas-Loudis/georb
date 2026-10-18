function [orbit_matrices] = orbit_matrices_init() 

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
 
 
%-------------------------------------------------------------------------- 
% Function: orbit_matrices 
%-------------------------------------------------------------------------- 
% Purpose: 
%  orbit_matrices initialisation  
%-------------------------------------------------------------------------- 
% Input arguments: 
% - orbit_matrices     : orbit_matrices strucuture array name 
% 
% Output arguments: 
% - orbit_matrices     : orbit_matrices strucuture with initial matrices 
% 
%-------------------------------------------------------------------------- 
% Dr. Thomas Loudis Papanikolaou                          28 september 2026 
%-------------------------------------------------------------------------- 
 
 
%-------------------------------------------------------------------------- 
% Orbits matrices to output structure array 
%-------------------------------------------------------------------------- 
mat_init = 0;

orbit_matrices.orbit_crf = mat_init; 
orbit_matrices.orbit_trf = mat_init; 
orbit_matrices.orbit_kepler = mat_init; 
 
orbit_matrices.state_transition_matrix = mat_init; 
orbit_matrices.sensitivity_matrix = mat_init; 
 
orbit_matrices.forces_acceleration = mat_init; 
orbit_matrices.gravity_gradient_matrix = mat_init; 
orbit_matrices.earth_orientation_matrix = mat_init; 
 
% orbit_matrices.observation_matrix = OBS_matrix; 
orbit_matrices.observation_residuals = mat_init; 
orbit_matrices.rms_observation_residuals = mat_init; 
 
% orbit_matrices.external_orbit_crf = extorb_crf; 
 
orbit_matrices.rms_orbital = mat_init; 
 
orbit_matrices.delta_orbital = mat_init; 
orbit_matrices.delta_orb_crf = mat_init; 
orbit_matrices.delta_kepler = mat_init; 
orbit_matrices.delta_orb_trf = mat_init;  
%-------------------------------------------------------------------------- 
