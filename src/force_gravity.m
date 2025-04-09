function [accel_vec, partials_r, partials_p, GM_Earth, radius_Earth, Cnm_tidefree, Snm_GFM] = force_gravity(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, legendre_functions_struct, gfm_struct_glob)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_gravity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Gravity Field' acceleration vector and partial derivatives
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'Earth_Gravity_Field';
[effect_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(effect_yn,'y');
if test == 1    
    % SHC Global array (structure)  
    struct_array = gfm_struct_glob;         
    % GM and radius
    GM_Earth   = struct_array.GM;
    radius_Earth = struct_array.radius;
    % Degree, Order truncation values
    if VEQ_mode_test == 1
    nm_values = struct_array.degree_veq;
    else
    nm_values = struct_array.degree_eqm;
    end
    degree_GFM = nm_values(1);
    order_GFM  = nm_values(2);
    % Spherical Harmonic coefficients matrices  
    Cnm_GFM = struct_array.Cnm;
    Snm_GFM = struct_array.Snm; 
else
    % Standards
    GM_Earth     = 3.9860044150e+14;
    radius_Earth = 6.3781363000e+06; 
    Cnm_GFM = 0;
    Snm_GFM = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Zero Tide" correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfm_tide_system = gfm_struct_glob.tide_system;
test = strcmp(gfm_tide_system,'zero_tide');   
if test == 1            
    % IERS Conventions 2010 : Solid Earth Tides Step3 : Correction of C20 for "Zero Tide" Gravity models
    %K20 = 0.29525;
    K20 = 0.30190;
    dC20_perm = (4.4228 * 10^(-8)) * (-0.31460) * K20;
    Cnm_tidefree = Cnm_GFM;
    Cnm_tidefree(2+1,0+1) = Cnm_GFM(2+1,0+1) - dC20_perm;
else
    Cnm_tidefree = Cnm_GFM;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'Gravity_field_terms';
[Gravity_field_terms_yn] = read_param_cfg(ORB_config,param_keyword);
test_central_grav = strcmp(Gravity_field_terms_yn,'central');
test_static_grav = strcmp(Gravity_field_terms_yn,'static');
test_timevar_grav = strcmp(Gravity_field_terms_yn,'time-variable');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central gravity field - Equivalent to the effect of the Coo term only 
if test_central_grav == 1
    [a_earth_x,a_earth_y,a_earth_z] = accel_gm(rGCRS,GM_Earth);
    % VEQ : Partial derivatives
    [Uearth] = pdv_acclgm(rGCRS,GM_Earth);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geopotential (full terms)
elseif test_static_grav == 1 || test_timevar_grav == 1   
    % Spherical Harmonic Synthesis start degree
    degree_min = 0;
    if VEQ_mode_test == 0

    % Acceleration in ITRS
    [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,degree_GFM,order_GFM,GM_Earth,radius_Earth,Cnm_tidefree,Snm_GFM, legendre_functions_struct, degree_min);
    ax = partials_xyz(1,1);
    ay = partials_xyz(2,1);
    az = partials_xyz(3,1);
    Uearth = zeros(3,3);
    elseif VEQ_mode_test == 1
        % Acceleration and Partials w.r.t. state vector in ITRS
        [partials_2nd_spher, partials_2nd_xyz, partials_1st_spher, partials_1st_xyz] = potential_partials_2nd(rITRS,degree_GFM,order_GFM,GM_Earth,radius_Earth,Cnm_tidefree,Snm_GFM, legendre_functions_struct, degree_min);        
        Uearth = partials_2nd_xyz;
        ax = partials_1st_xyz(1,1);
        ay = partials_1st_xyz(2,1);
        az = partials_1st_xyz(3,1);
    end
    % Transformation of acceleration from ITRS to the GCRS
    aGCRS = eopmatrix * [ax; ay; az];
    a_earth_x = aGCRS(1,1);
    a_earth_y = aGCRS(2,1);
    a_earth_z = aGCRS(3,1);
    %clear aGCRS ax ay az    
    if VEQ_mode_test == 1
        % 2nd order partials: Transformation from ITRS to GCRS
        % (df/dr)_Inertial = EOP(t) * (df/dr)_Terrestrial * inv( EOP(t) )
        Uearth = eopmatrix * Uearth * (eopmatrix)' ;       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    a_earth_x = 0;
    a_earth_y = 0;
    a_earth_z = 0;
    Uearth = zeros(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_earth =  [a_earth_x; a_earth_y; a_earth_z];
a_earth_vec = sqrt(a_earth(1,1)^2 + a_earth(2,1)^2 + a_earth(3,1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accel_vec = a_earth;
partials_r = Uearth;
partials_p = 0;
