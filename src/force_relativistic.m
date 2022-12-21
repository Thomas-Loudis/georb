function [accel_vec] = force_relativistic(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config, GM_Earth,GM_Sun,z_Sun)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_relativistic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Relativistic effects including Schwarzschild effect, Lense-Thirring
%  precession and geodesic effect or de Sitter precession
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


global relativistic_effects_glob 

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
% GM constant of Moon and Sun in m^3/sec^2 
GMsun = GM_Sun;
% Sun position vector in CRS
rgSun = z_Sun(1:3,1);
% Sun velocity vector in CRS
vgSun = z_Sun(4:6,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic corrections to satellite's acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct_array = relativistic_effects_glob;
Relativity_yn = struct_array.Relativity_yn;
Relativistic_effects_yn = struct_array.Relativistic_effects_yn;
ppn_beta = struct_array.PPN_beta;
ppn_gama = struct_array.PPN_gama;
cspeedlight = struct_array.speedoflight;

test = strcmp(Relativity_yn,'y');
if test == 1 && VEQ_mode_test == 0
    % State Vector
    zGCRS = Z_crs;
    % Earth state vector with respect to Sun
    Zearth = [-rgSun; -vgSun];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwarzschild effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_effect = strcmp(Relativistic_effects_yn(1),'y');
if test_effect == 1     
    effect_name = 'Schwarzschild';
    [a_Schwarzschild] = relativistic_effects(zGCRS,Zearth,GM_Earth,GMsun,cspeedlight,ppn_beta,ppn_gama,effect_name);
else
    a_Schwarzschild = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lense-Thirring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_effect = strcmp(Relativistic_effects_yn(2),'y');
if test_effect == 1     
    effect_name = 'Lense-Thirring';
    [a_LenseThirring] = relativistic_effects(zGCRS,Zearth,GM_Earth,GMsun,cspeedlight,ppn_beta,ppn_gama,effect_name);
else
    a_LenseThirring = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geodetic effect or de Sitter precession
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_effect = strcmp(Relativistic_effects_yn(3),'y');
if test_effect == 1     
    effect_name = 'deSitter';
    [a_deSitter] = relativistic_effects(zGCRS,Zearth,GM_Earth,GMsun,cspeedlight,ppn_beta,ppn_gama,effect_name);
else
    a_deSitter = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall relativistic effects corrections to acceleration vector
    a_relativistic = a_Schwarzschild + a_LenseThirring + a_deSitter;
    
else 
    a_relativistic = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accel_vec = a_relativistic;
