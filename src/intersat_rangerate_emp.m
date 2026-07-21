function [rangerate_emp_obs, PD_rangerate_emp_MATRIX] = intersat_rangerate_emp(mjd, mjd0, zA, zB, rangerate_param_matrix) 

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
% intersat_rangerate_emp:  Inter-satellite Range-Rate observations empirical parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Inter-satellite Range-Rate observations empirical corrections functional  
%  Empirical parameters include bias, bias drift and periodic terms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - filename:   KBR1b data file's name 
% 
% Output arguments: 
% - intersat_range:         Range matrix 
%   intersat_range = [MJD range] 
% - intersat_rangerate:     Range rate matrix 
%   intersat_rangerate = [MJD range_rate] 
% - intersat_rangeaccel:    Range acceleration data 
%   intersat_rangeaccel = [MJD range_acceleration] 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                           12 November 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
% GM: Earth gravitational constant in m^3/sec^2  
GM = 3.9860044150e+14; % 3986004.418*10^8; % IERS Conventions 2003    
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range-Rate Empirical parameters from aprior matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[Ncycles, d2] = size(rangerate_param_matrix); 
Ntime_arguments = 2; 
% Number of emprical parameters per period 
Nparam_EMP = d2 - Ntime_arguments; 
% Overall Number of emprical parameters for all periods 
Nparam_EMP_all = Nparam_EMP * Ncycles ; 
 
if Ncycles > 1 
intersat_rangerate_emp_step = rangerate_param_matrix(2,2) - rangerate_param_matrix(1,2);  
else 
intersat_rangerate_emp_step = inf;  
end 
 
% Partials Matrix for current epoch mjd initialisation to zeros 
PD_rangerate_emp_MATRIX = zeros(1,Nparam_EMP_all); 
 
t_mjd = mjd;  
t_sec = (mjd - fix(mjd)) * 86400; 
 
for i = 1 : Ncycles 
    mjd_cycle = rangerate_param_matrix(i,1); 
    sec_cycle = rangerate_param_matrix(i,2); 
    % if ( fix(mjd) == fix(mjd_cycle) ) && (t_sec >= sec_cycle) && (t_sec < sec_cycle + intersat_rangerate_emp_step)  
    if (mjd >= mjd_cycle) && (mjd < (mjd_cycle + intersat_rangerate_emp_step/86400) )  
        cycle_index = i; 
        rangerate_emp_param = rangerate_param_matrix(cycle_index, Ntime_arguments+1 : Ntime_arguments + Nparam_EMP)'; 
 
if Nparam_EMP == 6 
% Bias and Bias drift 
B1 = rangerate_emp_param(1,1); 
B2 = rangerate_emp_param(2,1); 
% Periodical cosine terms 
C1 = rangerate_emp_param(3,1); 
C2 = rangerate_emp_param(4,1); 
% Periodical sin terms 
S1 = rangerate_emp_param(5,1); 
S2 = rangerate_emp_param(6,1); 
 
elseif Nparam_EMP == 3 
% Bias and Bias drift 
B1 = rangerate_emp_param(1,1); 
B2 = 0; 
% Periodical cosine terms 
C1 = rangerate_emp_param(2,1); 
C2 = 0; 
% Periodical sin terms 
S1 = rangerate_emp_param(3,1); 
S2 = 0; 
end 
 
    end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
mjd0 = rangerate_param_matrix(1,1); 
% mjd0 = rangerate_param_matrix(cycle_index,1); 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Argument of Latitude referring to the midpoint between 2 orbits 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Time dependent argument of CPR cos and sin terms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Time difference since start of the orbit arc 
mjd_to = mjd0; 
delta_mjd = mjd - mjd_to; 
delta_tsec = delta_mjd * 24 * 3600; 
 
 
% Satellite Orbit 1 
rGCRS = zA(1:3,1); 
vGCRS = zA(4:6,1); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Keplerian elements 
[a_semiaxis,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_deg] = kepler(rGCRS,vGCRS,GM); 
u_rad = u_deg * (pi / 180); 
M_rad = M_deg * (pi / 180); 
 
% Orbital Period 
T_orb_period = 2*pi * sqrt(a_semiaxis^3 / GM); 
 
% Orbital Mean Motion 
n_motion = 2*pi / T_orb_period;   
 
% Time-variable argument of CPR terms 
cpr_time_argument = u_rad; 
cpr_time_argument = M_rad; 
cpr_time_argument = n_motion * delta_tsec; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
time_argument_sat1 = cpr_time_argument; 
 
 
% Satellite Orbit 2 
rGCRS = zB(1:3,1); 
vGCRS = zB(4:6,1); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Keplerian elements 
[a_semiaxis,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_deg] = kepler(rGCRS,vGCRS,GM); 
u_rad = u_deg * (pi / 180); 
M_rad = M_deg * (pi / 180); 
 
% Orbital Period 
T_orb_period = 2*pi * sqrt(a_semiaxis^3 / GM); 
 
% Orbital Mean Motion 
n_motion = 2*pi / T_orb_period;   
 
% Time-variable argument of CPR terms 
cpr_time_argument = u_rad; 
cpr_time_argument = M_rad; 
cpr_time_argument = n_motion * delta_tsec; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
time_argument_sat2 = cpr_time_argument; 
 
 
u_mid = (time_argument_sat1 + time_argument_sat2) / 2; 
t = delta_tsec; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Empirical Function as range-rate observation correction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rangerate_emp_obs = B1 + B2 * t + (C1 + C2*t) * cos(u_mid) + (S1 + S2*t) * sin(u_mid);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partials w.r.t. empirical parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PD_rangerate_emp = zeros(3,Nparam_EMP);  
PD_rangerate_emp = zeros(1,Nparam_EMP); 
 
if Nparam_EMP == 6 
 
% Bias and Bias drift 
PD_rangerate_emp(1,1) = 1; 
PD_rangerate_emp(1,2) = t; 
 
% Periodical cosine terms 
PD_rangerate_emp(1,3) = cos(u_mid); 
PD_rangerate_emp(1,4) = t * cos(u_mid); 
 
% Periodical sine terms 
PD_rangerate_emp(1,5) = sin(u_mid); 
PD_rangerate_emp(1,6) = t * sin(u_mid); 
 
elseif Nparam_EMP == 3 
 
% Bias and Bias drift 
PD_rangerate_emp(1,1) = 1; 
 
% Periodical cosine terms 
PD_rangerate_emp(1,2) = cos(u_mid); 
 
% Periodical sin terms 
PD_rangerate_emp(1,3) = sin(u_mid); 
 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partials Matrix for current epoch mjd 
PD_rangerate_emp_MATRIX(1, (cycle_index - 1) * Nparam_EMP + 1 : cycle_index * Nparam_EMP ) = PD_rangerate_emp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
