function [partials_p, accel_vec, partials_r] = force_gravity_estim(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, gfm_struct_glob, legendre_functions_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_gravity_estim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Gravity Field parameter estimation partial derivatives
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
% - partials_p:     Partials w.r.t. parameters 
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
grav_paramestim_yn = struct_array.param_estim_yn;
test_grav_paramestim_10 = strcmp(grav_paramestim_yn,'y');
if test_grav_paramestim_10 == 1
    % Gravity signal' C,S harmonics coefficients matrices to be estimated as corrections to the apriori gravity model Cnm, Snm
    delta_Cnm = struct_array.Cnm_estim; 
    delta_Snm = struct_array.Snm_estim; 

    % Gravity Field parameters to be estimated :: Maximum degree of coefficients
    Nrange          = struct_array.param_estim_degree;
    Nparam_grav_min = Nrange(1,1);
    Nparam_grav_max = Nrange(1,2);
    Nparam_grav     = struct_array.parameters_number;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field partials w.r.t. harmonics coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test_grav_paramestim_10 == 1    
    if VEQ_mode_test == 0    
        partials_p_coef = zeros(3,Nparam_grav);       
    elseif VEQ_mode_test == 1        
        % Partials w.r.t. gravity field parameters
        [partials_p_coef_TRS, partials_c_TRS, partials_s_TRS] = potential_partials_coef(rITRS,Nparam_grav_max,Nparam_grav_min,GM_Earth,radius_Earth,gfm_struct_glob);        
        % Transformation from Terrestrial (ITRS) to Celestial (ICRS) Reference frame
        partials_p_coef = eopmatrix * partials_p_coef_TRS;        
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    partials_p_coef = zeros(3,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partials_p = partials_p_coef;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration of the (time variable) gravity signal (estimable) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test_grav_paramestim_10 == 1    
    % Spherical Harmonic Synthesis start degree
    degree_min = Nparam_grav_min;
    % Degree and Order truncated values
    degree_trunc = Nparam_grav_max;
    order_trunc  = Nparam_grav_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
n_max = Nparam_grav_max;
% computation of spherical coordinates (in radians)
[lamda,phi,l] = lamda_phi(rITRS);
rdist = l;
% Normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,n_max);

% First-order derivatives of normalized associated Legendre functions
[dPnm_norm] = Legendre1ord(phi,n_max) ;

% Second-order derivatives of the Normalized Associated Legendre functions
[d2Pnm_norm] = Legendre2ord(phi,n_max) ;

legendre_functions_struct.Pnm_norm = Pnm_norm;
legendre_functions_struct.Pnm_norm_derivatives_1st = dPnm_norm;
legendre_functions_struct.Pnm_norm_derivatives_2nd = d2Pnm_norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    if VEQ_mode_test == 0    
        % Acceleration in ITRS
        [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,degree_trunc,order_trunc,GM_Earth,radius_Earth, delta_Cnm, delta_Snm, legendre_functions_struct, degree_min);
        ax = partials_xyz(1,1);
        ay = partials_xyz(2,1);
        az = partials_xyz(3,1);
        Uearth = zeros(3,3);
    elseif VEQ_mode_test == 1
        % Acceleration and Partials w.r.t. state vector in ITRS
        [partials_2nd_spher, partials_2nd_xyz, partials_1st_spher, partials_1st_xyz] = potential_partials_2nd(rITRS,degree_trunc,order_trunc,GM_Earth,radius_Earth, delta_Cnm, delta_Snm, legendre_functions_struct, degree_min);        
        Uearth = partials_2nd_xyz;
        ax = partials_1st_xyz(1,1);
        ay = partials_1st_xyz(2,1);
        az = partials_1st_xyz(3,1);
    end

    % Transformation of acceleration from ITRS to the GCRS
    aGCRS = eopmatrix * [ax; ay; az];
    a_grav_x = aGCRS(1,1);
    a_grav_y = aGCRS(2,1);
    a_grav_z = aGCRS(3,1);
    if VEQ_mode_test == 1
        % 2nd order partials: Transformation from ITRS to GCRS
        % (df/dr)_Inertial = EOP(t) * (df/dr)_Terrestrial * inv( EOP(t) )
        Uearth = eopmatrix * Uearth * (eopmatrix)' ;       
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    a_grav_x = 0;
    a_grav_y = 0;
    a_grav_z = 0;
    Uearth = zeros(3,3);
end
accel_grav_signal =  [a_grav_x; a_grav_y; a_grav_z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accel_vec = accel_grav_signal;
partials_r = Uearth;

