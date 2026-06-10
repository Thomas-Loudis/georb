function [rms_CS] = stat_harmonics_rms (CS_matrix, n_min, n_max, m_min) 

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
% Function: stat_harmonics_rms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Compute RMS of spherical harmonic coefficients numerical differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas Loudis Papanikolaou                                   28 June 2023 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Statistics :: RMS of delta coefficients 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
delta_CS = CS_matrix; 
N_delta = ( (n_max+1)^2 - (n_max+1) ) / 2 + (n_max+1); 
delta_matrix = zeros(N_delta,1); 
 
i_delta = 0; 
for n = n_min : n_max 
    for m = m_min : n 
        i_delta = i_delta + 1; 
        delta_matrix(i_delta,1) = delta_CS(n+1,m+1); 
    end 
end 
% i_delta 
delta_matrix_out = delta_matrix(1:i_delta,1); 
 
% RMS of delta coefficients 
rms_CS = rms(delta_matrix_out); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
