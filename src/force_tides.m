function [accel_tides, accel_atm_tides, partials_r] = force_tides(mjd,Z_crs,Rtrs2crs, EQ_mode,ORB_config, GM_Earth,radius_Earth,GM_Moon,r_Moon,GM_Sun,r_Sun, eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_tides 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Earth Tides effects including:
%  - Solid Earth Tides 
%  - Ocean Tides 
%  - Solid Earth Pole Tide 
%  - Ocean Pole Tide
%  - Atmospheric Tides
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


global ocean_tides_struct_glob 
global aod_tides_struct_glb

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
GMmoon = GM_Moon;
GMsun = GM_Sun;
% Body-Fixed positions of Sun and Moon (Transformation : GCRS to ITRS)
rgMoon_ITRS = (eopmatrix)' * r_Moon;
rgSun_ITRS = (eopmatrix)' * r_Sun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'Tides_effects';
[Tides_effects] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(Tides_effects,'y');
if test == 1    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-independent : Step 1 (IERS Conventions 2010)
    param_keyword = 'solid_earth_tides_1_non_freq';
    [solid_earth_tides_1_non_freq] = read_param_cfg(ORB_config,param_keyword);
    test = strcmp(solid_earth_tides_1_non_freq,'y');
    if test == 1
        [dCnm_solid1,dSnm_solid1] = tides_solid1(GM_Earth,radius_Earth,GMmoon,rgMoon_ITRS,GMsun,rgSun_ITRS);
         % Solid Earth Tides (1) acceleration vector (in ITRS)
        %[ax,ay,az] = accel_tide(rITRS,4,4,GM_Earth,radius_Earth,dCnm_solid1,dSnm_solid1);
        [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,4,4,GM_Earth,radius_Earth,dCnm_solid1,dSnm_solid1);
        ax = partials_xyz(1,1);
        ay = partials_xyz(2,1);
        az = partials_xyz(3,1);        
        % Transformation of acceleration from ITRS to the GCRS
        a_solid1 = eopmatrix * [ax; ay; az];
    else
        dCnm_solid1 = zeros(4+1,4+1);
        dSnm_solid1 = zeros(4+1,4+1); 
        a_solid1 = [0 0 0]';
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-dependent : Step 2 (IERS Conventions 2010)
    param_keyword = 'solid_earth_tides_2_freq';
    [solid_earth_tides_2_freq] = read_param_cfg(ORB_config,param_keyword);
    test = strcmp(solid_earth_tides_2_freq,'y');
    if test == 1
        [dCnm_solid2,dSnm_solid2] = tides_solid2(mjd,eop,dpint);
         % Solid Earth Tides (2) acceleration vector (in ITRS)
        %[ax,ay,az] = accel_tide(rITRS,2,2,GM_Earth,radius_Earth,dCnm_solid2,dSnm_solid2);
        [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,2,2,GM_Earth,radius_Earth,dCnm_solid2,dSnm_solid2);
        ax = partials_xyz(1,1);
        ay = partials_xyz(2,1);
        az = partials_xyz(3,1);
        % Transformation of acceleration from ITRS to the GCRS
        a_solid2 = eopmatrix * [ax; ay; az];
    else
         dCnm_solid2 = zeros(2+1,2+1);
         dSnm_solid2 = zeros(2+1,2+1);
        a_solid2 = [0 0 0]';
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param_keyword = 'ocean_tides';
    [ocean_tides_yn] = read_param_cfg(ORB_config,param_keyword);
    test = strcmp(ocean_tides_yn,'y');
    if test == 1
        % Global array (structure) of Atmospheric Tides   
        struct_array = ocean_tides_struct_glob;         
        % Spherical Harmonic coefficients matrices  
        if VEQ_mode_test == 1
        nm_values = struct_array.degree_veq;
        else
        nm_values = struct_array.degree_eqm;
        end
        degree_trunc = nm_values(1);
        order_trunc  = nm_values(2);
        %
        ocean_tides_delaunay_doodson_multipliers = struct_array.delaunay_doodson_multipl;
        ocean_tides_dCnm_plus  = struct_array.dCnm_plus;
        ocean_tides_dSnm_plus  = struct_array.dSnm_plus;
        ocean_tides_dCnm_minus = struct_array.dCnm_minus;
        ocean_tides_dSnm_minus = struct_array.dSnm_minus;         
        % Ocean Tide model
        [dCnm_ocean,dSnm_ocean] = tides_ocean(degree_trunc,order_trunc,mjd,eop,dpint,ocean_tides_delaunay_doodson_multipliers,ocean_tides_dCnm_plus,ocean_tides_dSnm_plus,ocean_tides_dCnm_minus,ocean_tides_dSnm_minus);
        % Ocean Tides Acceleration vector (in ITRS)
        [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,degree_trunc,order_trunc,GM_Earth,radius_Earth,dCnm_ocean,dSnm_ocean);
        ax = partials_xyz(1,1);
        ay = partials_xyz(2,1);
        az = partials_xyz(3,1);        
        % Transformation of acceleration from ITRS to the GCRS
        a_ocean_GCRS = eopmatrix * [ax; ay; az];

    else
        dCnm_ocean = zeros(10+1,10+1);
        dSnm_ocean = zeros(10+1,10+1);
        a_ocean_GCRS = [0 0 0]';
    end
    % Ocean Tides acceleration vector (in celestial reference frame GCRS)
    a_ocean = a_ocean_GCRS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pole Tide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Pole Tide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'solid_earth_pole_tide';
[solid_earth_pole_tide_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(solid_earth_pole_tide_yn,'y');
if test == 1
[dCnm_se_pole_tide, dSnm_se_pole_tide] = tides_pole_solidearth(mjd,eop,dpint);
% Acceleration vector (in ITRS)
%[ax,ay,az] = accel_tide(rITRS,2,2,GM_Earth,radius_Earth,dCnm_se_pole_tide,dSnm_se_pole_tide);
[partials_rpl, partials_xyz] = potential_partials_1st(rITRS,2,2,GM_Earth,radius_Earth,dCnm_se_pole_tide,dSnm_se_pole_tide);
ax = partials_xyz(1,1);
ay = partials_xyz(2,1);
az = partials_xyz(3,1);
% Transformation of acceleration from ITRS to the GCRS
a_se_pole_tide = eopmatrix * [ax; ay; az];
else
a_se_pole_tide = [0 0 0]';
dCnm_se_pole_tide = zeros(3,3);
dSnm_se_pole_tide = zeros(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Pole Tide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'ocean_pole_tide';
[ocean_pole_tide_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(ocean_pole_tide_yn,'y');
if test == 1
[dCnm_ocean_pole_tide, dSnm_ocean_pole_tide] = tides_pole_ocean(mjd,eop,dpint);
% Acceleration vector (in ITRS)
[partials_rpl, partials_xyz] = potential_partials_1st(rITRS,2,2,GM_Earth,radius_Earth,dCnm_ocean_pole_tide,dSnm_ocean_pole_tide);
ax = partials_xyz(1,1);
ay = partials_xyz(2,1);
az = partials_xyz(3,1);
% Transformation of acceleration from ITRS to the GCRS
a_ocean_pole_tide = eopmatrix * [ax; ay; az];
else
a_ocean_pole_tide = [0 0 0]';
dCnm_ocean_pole_tide = zeros(3,3);
dSnm_ocean_pole_tide = zeros(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tides SHC array :: corrections to Cnm,Snm   
    % Solid Earth Tides corrections array
    [dCnm_solid,dSnm_solid] = tides_add2(dCnm_solid1,dSnm_solid1,dCnm_solid2,dSnm_solid2,-1);
    % Tides corrections array including Ocean Tides
    [dCnm_tide,dSnm_tide] = tides_add2(dCnm_ocean,dSnm_ocean,dCnm_solid,dSnm_solid,-1);
    %[dCnm_tide,dSnm_tide] = tides_add2(dCnm_solid,dSnm_solid,dCnm_ocean,dSnm_ocean,-1);    
    % Include the Pole Tide corrections
    dCnm_tide(2+1,1+1) = dCnm_tide(2+1,1+1) + dCnm_se_pole_tide(2+1,1+1) + dCnm_ocean_pole_tide(2+1,1+1);
    dSnm_tide(2+1,1+1) = dSnm_tide(2+1,1+1) + dSnm_se_pole_tide(2+1,1+1) + dSnm_ocean_pole_tide(2+1,1+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Tides acceleration summary
a_tides = a_solid1 + a_solid2 + a_ocean + a_se_pole_tide + a_ocean_pole_tide;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    dCnm_tide = 0;
    dSnm_tide = 0;
    a_tides = [0 0 0]';
    Utides_icrs = zeros(3,3);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'atm_tides';
[atm_tides_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(atm_tides_yn,'y');
if test == 1
    % Global array (structure) of Atmospheric Tides   
    struct_array = aod_tides_struct_glb;         
    % Atmospheric Tides Spherical Harmonic coefficients matrices  
    atmtides_GM      = struct_array.GM;
    atmtides_radius  = struct_array.radius;
    %atm_tides_degree = struct_array.degree
    if VEQ_mode_test == 1
    nm_values = struct_array.degree_veq;
    else
    nm_values = struct_array.degree_eqm;
    end
    degree_trunc = nm_values(1);
    order_trunc  = nm_values(2);
    %
    atmtides_delaunay_doodson_multipliers = struct_array.delaunay_doodson_multipl;
    atm_tides_dCnm_plus  = struct_array.dCnm_plus;
    atm_tides_dSnm_plus  = struct_array.dSnm_plus;
    atm_tides_dCnm_minus = struct_array.dCnm_minus;
    atm_tides_dSnm_minus = struct_array.dSnm_minus;
    % Matrices dCnm_atm_tides, dSnm_atm_tides
    [dCnm_atm_tides,dSnm_atm_tides] = tides_ocean(degree_trunc,order_trunc,mjd,eop,dpint,atmtides_delaunay_doodson_multipliers,atm_tides_dCnm_plus,atm_tides_dSnm_plus,atm_tides_dCnm_minus,atm_tides_dSnm_minus); 
    % Atmospheric Tides Acceleration vector (in ITRS)
    %[ax,ay,az] = accel_tide(rITRS,degree_trunc,order_trunc,atmtides_GM,atmtides_radius,dCnm_atm_tides,dSnm_atm_tides);
    [partials_rpl, partials_xyz] = potential_partials_1st(rITRS,degree_trunc,order_trunc,atmtides_GM,atmtides_radius,dCnm_atm_tides,dSnm_atm_tides);
    ax = partials_xyz(1,1);
    ay = partials_xyz(2,1);
    az = partials_xyz(3,1);
    % Transformation of acceleration from ITRS to the GCRS
    a_atm_tides_crf = eopmatrix * [ax; ay; az]; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tides corrections array including Atmospheric Tides
    [dCnm_tides_sum,dSnm_tides_sum] = tides_add2(dCnm_atm_tides,dSnm_atm_tides,dCnm_tide,dSnm_tide,-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
else
    a_atm_tides_crf = [0; 0; 0];
    dCnm_tides_sum = dCnm_tide;
    dSnm_tides_sum = dSnm_tide;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ : Partial derivatives of tides acceleration w.r.t. state vector in ITRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'Tides_effects';
[Tides_effects] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(Tides_effects,'y');
if test == 1       
Utides_icrs = zeros(3,3);    
if VEQ_mode_test == 1
    [d1 d2] = size(dCnm_tides_sum);
    n_tides = d1 - 1;
    m_tides = d2 - 1;
    % Geopotential 2nd order partial derivatives 
    [partials_2nd_spher, partials_2nd_xyz, partials_1st_spher, partials_1st_xyz] = potential_partials_2nd(rITRS,n_tides,m_tides,GM_Earth,radius_Earth,dCnm_tides_sum,dSnm_tides_sum);
    Utides_itrs = partials_2nd_xyz;
    % Transformation from ITRS to GCRS
    % (df/dr)_Inertial = EOP(t) * (df/dr)_Terrestrial * inv( EOP(t) )
    Utides_icrs = eopmatrix * Utides_itrs * (eopmatrix)' ;
else
    Utides_icrs = zeros(3,3);    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accel_tides = a_tides;
accel_atm_tides = a_atm_tides_crf;
partials_r = Utides_icrs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

