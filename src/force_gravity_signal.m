function [accel_vec] = force_gravity_signal(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, gfm_struct_glob, legendre_functions_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_gravity_estim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Gravity (temporal) signal acceleration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
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
% - accel_vec:      Acceleration vector cartesian components in Celestial reference frame (GCRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Thomas Loudis Papanikolaou                                     
% Created: 21 March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% Gravitational effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field structure array  
struct_array = gfm_struct_glob;         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'Earth_Gravity_Field';
[effect_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(effect_yn,'y');
if test == 1    
    % GM and radius
    GM_Earth   = struct_array.GM;
    radius_Earth = struct_array.radius;    
else
    % Standards
    GM_Earth     = 3.9860044150e+14;
    radius_Earth = 6.3781363000e+06; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters estimation y/n 
gravity_signal_yn = struct_array.gravity_signal_yn;
test_effect_01 = strcmp(gravity_signal_yn,'y');
if test_effect_01 == 1
    % Gravity signal' C,S harmonics coefficients matrices to be estimated as corrections to the apriori gravity model Cnm, Snm
    delta_Cnm = struct_array.Cnm_signal; 
    delta_Snm = struct_array.Snm_signal; 
    
    % Gravity Field parameters to be estimated :: Maximum degree of coefficients
    Nrange          = struct_array.param_estim_degree;
    Nparam_grav_min = Nrange(1,1);
    Nparam_grav_max = Nrange(1,2);
    Nparam_grav     = struct_array.parameters_number;

    degree_order_signal = struct_array.degree_signal; 
    degree_signal = degree_order_signal(1,1);
    order_signal  = degree_order_signal(1,2);
    
    % Spherical Harmonic Synthesis start degree
    degree_min = Nparam_grav_min;
    % Degree and Order truncated values
    % degree_trunc = Nparam_grav_max;
    % order_trunc  = Nparam_grav_max;
    degree_trunc = degree_signal;
    order_trunc  = order_signal;
    
    % Degree, Order truncation values
    if VEQ_mode_test == 1
    nm_values = struct_array.degree_veq;

    degree_GFM = nm_values(1);
    order_GFM  = nm_values(2);    
    if degree_trunc > degree_GFM
        degree_trunc = degree_GFM;
        order_trunc  = order_GFM;
    end    
    % else
    % nm_values = struct_array.degree_eqm;
    end
    % degree_GFM = nm_values(1);
    % order_GFM  = nm_values(2);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration of the (time variable) gravity signal (estimable) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acceleration in ITRS
    [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,degree_trunc,order_trunc,GM_Earth,radius_Earth, delta_Cnm, delta_Snm, legendre_functions_struct, degree_min);
    ax = partials_xyz(1,1);
    ay = partials_xyz(2,1);
    az = partials_xyz(3,1);

    % Transformation of acceleration from ITRS to the GCRS
    aGCRS = eopmatrix * [ax; ay; az];
    a_grav_x = aGCRS(1,1);
    a_grav_y = aGCRS(2,1);
    a_grav_z = aGCRS(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    a_grav_x = 0;
    a_grav_y = 0;
    a_grav_z = 0;
    % Uearth = zeros(3,3);
end
accel_grav_signal =  [a_grav_x; a_grav_y; a_grav_z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accel_vec = accel_grav_signal;
% partials_r = Uearth;

