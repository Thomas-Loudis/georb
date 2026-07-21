function [rangerate_empirical_struct] = prm_intersat_empirical (orbit_model_struct) 

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
% Function: prm_intersat_empirical 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbit_model_struct          : Orbit model structure array  
% 
% Output arguments: 
% - rangerate_empirical_struct  : Inter-satellite range-rate empirical parameters structure array  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                           19 November 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range-Rate Empirical Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rangerate_empirical_function_01 = orbit_model_struct.intersat_rangerate_empirical_01; 
 
% rangerate_empirical_struct.rangerate_empirical_function_01 = 1 
rangerate_empirical_struct.rangerate_empirical_function_01 = rangerate_empirical_function_01; 
 
% Step interval in seconds 
rangerate_empirical_struct.rangerate_empirical_step = 10800; 
 
% Observation function model: Number of empirical parameters per cycle e.g  
% bias and bias dirft and periodical terms (sine and cosine) coefficients 
% and their derivatives 
% rangerate_empirical_struct.rangerate_empirical_param_number = 6 
rangerate_empirical_struct.Nparam_rangerate_cycle  = 6; 
 
orbit_arc_length = orbit_model_struct.orbit_arc_length_sec;  
IC_MJDo          = orbit_model_struct.IC_MJD;  
 
% Parameters Matrix formulation, initialisation, preallocation 
intersat_rangerate_emp_step = rangerate_empirical_struct.rangerate_empirical_step;  
Nparam_obs_rangerate = rangerate_empirical_struct.Nparam_rangerate_cycle; 
 
 
Ncycles = fix(orbit_arc_length/intersat_rangerate_emp_step) + 1;  % % -1 - (ceil(pulse_offset/pulse_step))  % fix(orbit_arc_length/pulse_step)+1 - (ceil(pulse_offset/pulse_step)); % ceil(orbit_arc_length/pulse_step) - (ceil(pulse_offset/pulse_step))   
 
 
rangerate_empirical_parameters_matrix = zeros(Ncycles, 2 + Nparam_obs_rangerate); 
mjd0 = IC_MJDo; 
sec00 = (mjd0 - fix(mjd0)) * 86400; 
for i = 1 : Ncycles 
rangerate_empirical_parameters_matrix(i,1) = mjd0 + (i-1) * intersat_rangerate_emp_step / 86400; 
rangerate_empirical_parameters_matrix(i,2) = sec00+ (i-1) * intersat_rangerate_emp_step; 
end 
rangerate_empirical_parameters_apriori = rangerate_empirical_parameters_matrix; 
 
 
rangerate_empirical_struct.rangerate_empirical_parameters  = rangerate_empirical_parameters_apriori; 
% rangerate_empirical_struct.Nparam_rangerate_cycle  = Nparam_obs_rangerate; 
rangerate_empirical_struct.Nparam_rangerate = Nparam_obs_rangerate * Ncycles; 
rangerate_empirical_struct.rangerate_empirical_cycles_number  = Ncycles; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
