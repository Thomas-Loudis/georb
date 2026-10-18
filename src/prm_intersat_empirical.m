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
 
 
%-------------------------------------------------------------------------- 
% Range-Rate Empirical Parameters 
%-------------------------------------------------------------------------- 
rangerate_empirical_struct = orbit_model_struct.rangerate_empirical_struct ; 

intersat_rangerate_empirical_yn     = rangerate_empirical_struct.intersat_rangerate_empirical_yn;
% intersat_rangerate_empirical_step   = rangerate_empirical_struct.intersat_rangerate_empirical_step;

test_status_yn = strcmp(intersat_rangerate_empirical_yn,'y');    
if test_status_yn == 1    
rangerate_empirical_function_01 = 1; 
else 
rangerate_empirical_function_01 = 0;
end     

% rangerate_empirical_function_01 = orbit_model_struct.intersat_rangerate_empirical_01;  
rangerate_empirical_struct.rangerate_empirical_function_01 = rangerate_empirical_function_01; 
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Observation function model: Number of empirical parameters per cycle   
% - Set Nparam_rangerate_cycle = 6 :: bias and bias drift and
%                                     periodical (sine and cosine) terms coefficients and their derivatives 
% - Set Nparam_rangerate_cycle = 3 :: bias and periodical (sine and cosine) terms coefficients
%                                     (Drift and derivatives are not estimated) 
% - Set Nparam_rangerate_cycle = 2 :: periodical (sine and cosine) terms coefficients
%                                     (derivatives are not estimated) 
rangerate_empirical_struct.Nparam_rangerate_cycle  = 6;
%-------------------------------------------------------------------------- 
 
orbit_arc_length = orbit_model_struct.orbit_arc_length_sec;  
IC_MJDo          = orbit_model_struct.IC_MJD;  

%-------------------------------------------------------------------------- 
% Step interval in seconds 
%-------------------------------------------------------------------------- 
intersat_rangerate_empirical_step   = rangerate_empirical_struct.intersat_rangerate_empirical_step;

% rangerate_empirical_struct.rangerate_empirical_step = 10800; 
rangerate_empirical_struct.rangerate_empirical_step = intersat_rangerate_empirical_step; 

intersat_rangerate_emp_step = rangerate_empirical_struct.rangerate_empirical_step;  
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Parameters Matrix formulation, initialisation, preallocation 
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
%-------------------------------------------------------------------------- 
  
rangerate_empirical_struct.rangerate_empirical_parameters  = rangerate_empirical_parameters_apriori; 
% rangerate_empirical_struct.Nparam_rangerate_cycle  = Nparam_obs_rangerate; 
rangerate_empirical_struct.Nparam_rangerate = Nparam_obs_rangerate * Ncycles; 
rangerate_empirical_struct.rangerate_empirical_cycles_number  = Ncycles; 
