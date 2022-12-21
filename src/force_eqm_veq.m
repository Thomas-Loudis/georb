function [ax,ay,az,pdv_acc,pdv_acc_param] = force_eqm_veq(z,eop,dpint, EQ_mode)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force model: acceleration vector and partials
% Function: force_eqm_veq 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the acceleration vector and force' partial derivatives
%  w.r.t. state vector and parametrs to be estimated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - z:       Epoch's MJD and State vector
%            z = [mjd r' v']            
%            mjd: Modified Julian Day number (including fraction of the day) in TT (Terrestrial Time)
%            r: Position vector in Celestial Reference System (GCRS)
%            v: Velocity vector in Celestial Reference System (GCRS)
% - eop:     Earth Orientation Parameters (EOP) data that are required for
%            the orbit arc length
% - dpint:   Number of data points (days) that are required for the EOP
%            interpolation to the computation epoch
% - EQ_mode: Inegral Equations mode: Equation of Motion or Variational Equations 
%
% Output arguments:
% - ax,ay,az        : Acceleration vector' cartesian components in GCRS
% - pdv_acc         : Force' partial derivatives w.r.t. state vector 
% - pdv_acc_param   : Force' partial derivatives w.r.t. parameters of interest 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                      November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 01/05/2011  upgrade for including planetary and tides perturbations
% 23/06/2012  New approach to Ocean Tides perturbations with FES2004 model
%             for minimizing computation time. 
%             Use of upgraded functions tides_ocean2.m & tides_fes2004_2.m
% 26/10/2022  Thomas Loudis Papanikolaou
%             Code upgrade and merge of functions accel.m and veq_accl.m 
%             into one function as force_eqm_veq.m
% 15/12/2022  Thomas Loudis Papanikolaou
%             Code upgrade; Writing functions force_xxxxx.m per effect
%             to be called by force_eqm_veq.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global ORB_config_struct_glob
global Nmodel_PARAM_ESTIM_glob Nparam_GLOB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations mode: EQM or VEQ
VEQ_mode_test = strcmp(EQ_mode,'VEQ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJD of computation epoch
mjd = z(1,1);
% JD of computation epoch
jd = mjd + 2400000.5;
% Position vector in GCRS
rGCRS = z(1,2:4)';
% Velocity vector in GCRS
vGCRS = z(1,5:7)';
% State Vector in GCRS
zGCRS = [rGCRS; vGCRS];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Celestial (GCRS) to Terrestrial (ITRS) transformation
[eopmatrix,deopmatrix] = trs2crs(mjd,eop,dpint);
% Position vector in ITRS
rITRS = (eopmatrix)' * rGCRS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master Orbit configuration modelling matrix
ORB_config = ORB_config_struct_glob;
% State Vector
Z_crs = [rGCRS; vGCRS];
% Terrestrial to Celstial refrence frame transformation matrix
Rtrs2crs = eopmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accel_vec, partials_r, partials_p, GM_Earth, radius_Earth, Cnm_tidefree, Snm_GFM] = force_gravity(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config);
a_earth = accel_vec;
Uearth = partials_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sun, Moon and Planets orbital perturbations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C20 = Cnm_tidefree(2+1,0+1);
[accel_bodies, accel_indirectJ2, partials_r, Moon_Z_crs, Sun_Z_crs, GM_Moon,GM_Sun] = force_planets(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, GM_Earth, radius_Earth, C20);
a_bodies = accel_bodies;
a_indirectJ2 = accel_indirectJ2;
pdv_3rdbodies = partials_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tides: 
% Solid Earth Tides (Frequency-independent & Frequency-dependent terms)
% Ocean Tides
% Solid Earth Pole Tide 
% Ocean Pole Tide 
% Atmospheric Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_Moon = Moon_Z_crs(1:3,1);
r_Sun  = Sun_Z_crs(1:3,1);
[accel_tides, accel_atm_tides, partials_r] = force_tides(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config, GM_Earth,radius_Earth,GM_Moon,r_Moon,GM_Sun,r_Sun, eop,dpint);
a_tides = accel_tides;
a_atm_tides_crf = accel_atm_tides;
Utides_icrs = partials_r;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic Effects as corrections to satellite's acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accel_vec] = force_relativistic(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config, GM_Earth,GM_Sun,Sun_Z_crs);
a_relativistic = accel_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmosphere and Ocean De-Aliasing (AOD) effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accel_vec] = force_aod(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config);
a_AOD_crf = accel_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall gravitational effects components
a_grav = a_earth + a_bodies + a_indirectJ2 + a_tides + a_relativistic + a_AOD_crf + a_atm_tides_crf;
pdv_grav = Uearth + Utides_icrs + pdv_3rdbodies;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Gravitational effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'non_gravitational_forces';
[non_gravitational_forces] = read_param_cfg(ORB_config_struct_glob,param_keyword);
test_nongrav = strcmp(non_gravitational_forces,'y');
if test_nongrav == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Accelerometer data processing and calibration modelling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
[accel_vec, partials_r, partials_p] = force_acc(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config);
facc_ICRF = accel_vec;
pdv_nongrav = partials_r;
PD_ACC_Cal_Param = partials_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

elseif test_nongrav == 0
    facc_ICRF = [0; 0; 0];
    % VEQ pdv
    pdv_nongrav = zeros(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall non-gravitational effects components
a_nongrav = facc_ICRF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical forces modelling based on Cycle-Per-Revolution acceleration terms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accel_vec, partials_r, partials_p] = force_cpr(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, GM_Earth);
a_emp = accel_vec;
PD_emp_Z = partials_r;
PD_accemp_P = partials_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Accelerations :: Pulses or Piecewise accelerations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accel_vec, partials_r, partials_p] = force_empaccel(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config);
f_pulses = accel_vec;
PD_pulses_param = partials_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summation acceleration vector acting on the satellite/object in orbit
% a_sum = a_grav + a_nongrav;
a_sum = a_grav + a_nongrav + a_emp + f_pulses;
ax = a_sum(1,1);
ay = a_sum(2,1);
az = a_sum(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives matrix of overall acceleration w.r.t. time 
pdv_acc = pdv_grav + pdv_nongrav + PD_emp_Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives matrix of overall acceleration w.r.t. unknown parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdv_acc_param = zeros(3,Nparam_GLOB);
if VEQ_mode_test == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Nparam_GLOB > 0)
Nmodel_PARAM_ESTIM = Nmodel_PARAM_ESTIM_glob;
Nparam = 0;    

% Empirical Forces :: CPR terms 
if Nmodel_PARAM_ESTIM(1) == 1
    [n1 n2] = size(PD_accemp_P);
    for i_param = 1 : n2
        Nparam = Nparam + 1;
        pdv_acc_param(:,Nparam) = PD_accemp_P(:,i_param);
    end
end

% Accelerometer calibration parameters
if Nmodel_PARAM_ESTIM(2) == 1
    [n1 n2] = size(PD_ACC_Cal_Param);
    for i_param = 1 : n2
        Nparam = Nparam + 1;
        pdv_acc_param(:,Nparam) = PD_ACC_Cal_Param(:,i_param);
    end
end

% Empirical Accelerations :: Piecewise acclerations or Pulses 
if Nmodel_PARAM_ESTIM(3) == 1
    [n1 n2] = size(PD_pulses_param);
    for i_param = 1 : n2
        Nparam = Nparam + 1;
        pdv_acc_param(:,Nparam) = PD_pulses_param(:,i_param);
    end
end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
