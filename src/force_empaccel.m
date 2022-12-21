function [accel_vec, partials_r, partials_p] = force_empaccel(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_empaccel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Empirical accelerations modelling based on piecewise constant
%  accelerations or pulses (instant velocity changes) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - mjd:            Epoch's Modified Julian Day number including fraction of the day) in TT (Terrestrial Time)
% - Z_crs:          State vector in Celestial Reference Frame (GCRS)   
%                   z = [r' v']
%                   r:   Position Vector (m) 
%                   v:   Velocity vector (m/sec) 
% - Rtrs2crs:       Tranformation matrix: Terrestrial Reference Frame to Celestial Reference Frame 
% - EQ_mode:        Inegral Equations mode: Equation of Motion or Variational Equations 
% - Configuration array:    Orbit master configuration structure array
%
% Output arguments:
% - accel_vec:      Acceleration vector cartesian components in GCRS
% - partials_r:     Partials w.r.t. position vector 
% - partials_p:     Partials w.r.t. parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Thomas Loudis Papanikolaou                                     
% Created: 15 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Code extracted from function force_eqm_veq. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global pulses_stoch_accel_glob 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations mode: EQM or VEQ
VEQ_mode_test = strcmp(EQ_mode,'VEQ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JD of computation epoch
jd = mjd + 2400000.5;
% Position vector in GCRS
rGCRS = Z_crs(1:3,1);
% Velocity vector in GCRS
vGCRS = Z_crs(4:6,1); 
% Terrestrial (ITRS) to Celestial (GCRS) 
eopmatrix = Rtrs2crs;
% Position vector in ITRS
rITRS = (eopmatrix)' * rGCRS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Accelerations :: Pulses or Piecewise accelerations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulses_accel_struct = pulses_stoch_accel_glob;
PULSES_estim_yn = pulses_accel_struct.effect_01;
test = strcmp(PULSES_estim_yn,'y');
if test == 1
    Sec_00 = (mjd - fix(mjd)) * 86400;    
    pulses_matrix      = pulses_accel_struct.acceleration_matrix;
    pulses_axes_vec_01 = pulses_accel_struct.parameters_01;
    N_param_pulses     = pulses_accel_struct.parameters_number;
    accel_duration     = pulses_accel_struct.duration;
    emp_accel_type = pulses_accel_struct.parameters_type;
    [f_pulses, PD_pulses_param, PD_pulses_r, PD_pulses_v] = stoch_accel(rGCRS, vGCRS, mjd, Sec_00, accel_duration, pulses_matrix, pulses_axes_vec_01, N_param_pulses,emp_accel_type);
else    
    f_pulses = [0 0 0]';
    PD_pulses_r = zeros(3,3);
    PD_pulses_param = zeros(3,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments
accel_vec  = f_pulses;
partials_r = PD_pulses_r;
partials_p = PD_pulses_param;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
