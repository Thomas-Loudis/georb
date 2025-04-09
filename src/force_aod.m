function [accel_vec] = force_aod(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config, legendre_functions_struct, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_aod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
% Atmosphere and Ocean De-Aliasing (AOD) effects' acceleration
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
% - accel_bodies:       Acceleration vector cartesian components in GCRS; Sum vector of all 3rd bodies perturbation
% - accel_indirectJ2:   Indirect J2 effect perturbation 
% - partials_r:         Partials w.r.t. position vector 
% - Moon_Z_crs:         Moon state vector in Celestial reference frame 
% - Sun_Z_crs:          Sun state vector in Celestial reference frame 
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
% Atmosphere and Ocean De-Aliasing (AOD) effects' acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'aod_effects';
[aod_effects_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(aod_effects_yn,'y');
if test == 1  &&  VEQ_mode_test == 0  
    % Forces model' structure array 
    aod_struct_glob = orbit_model_struct.aod_effects;
    struct_array = aod_struct_glob;         
    % Spherical Harmonic coefficients matrices  
    GM_shc      = struct_array.GM;
    radius_shc  = struct_array.radius;
    if VEQ_mode_test == 1
    nm_values = struct_array.degree_veq;
    else
    nm_values = struct_array.degree_eqm;
    end
    degree_trunc = nm_values(1);
    order_trunc  = nm_values(2);
    % SHC matrices    
    Cnm_shc  = struct_array.Cnm;
    Snm_shc  = struct_array.Snm;
    % Acceleration vector computation based on sherical harmonics series
    [ax_AOD,ay_AOD,az_AOD] = accel_aod(rITRS,degree_trunc,order_trunc,GM_shc,radius_shc, Cnm_shc, Snm_shc, mjd, legendre_functions_struct, struct_array);
    % Transformation of acceleration from ITRS to the GCRS
    a_AOD_crf = eopmatrix * [ax_AOD; ay_AOD; az_AOD]; 
else
    a_AOD_crf = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accel_vec = a_AOD_crf;
