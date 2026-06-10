function [orbit_model_matrix] = orbit_model_gravity (grav_field_paramestim_yn,orbit_model_matrix) 

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
% Function: orbit_model_gravity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Parameters/Variables for gravity model and gravity field determination  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                               7 April 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 30/04/2025, Thomas Loudis Papanikolaou 
%             Code minor modifications  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Earth Gravity Field 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if 1 == 0 
% Gravity model data file (spherical harmonic coefficients) 
[GM,ae,Cnm,Snm,sCnm,sSnm, n_max_eqm, m_max_eqm, n_max_veq, m_max_veq, tide_system, gravity_field_struct] = prm_grav_model(orbit_config_matrix); 
N_param_GRAV = gravity_field_struct.parameters_number; 
orbit_model_matrix.gravity_field = gravity_field_struct; 
orbit_model_matrix.GM_Earth = GM; 
end 
 
% Gravity Field structure matrix 
gravity_field_struct = orbit_model_matrix.gravity_field; 
% Update Gravity Field parameters estimation 
gravity_field_struct.param_estim_yn = grav_field_paramestim_yn; 
 
% Update orbit_model_matrix 
orbit_model_matrix.gravity_field = gravity_field_struct; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Update orbit_model_struct :: Forces effects and Number of parameters to be estimated  
[orbit_model_matrix] = orbit_model_force_parameters (orbit_model_matrix); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
