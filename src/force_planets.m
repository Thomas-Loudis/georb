function [accel_bodies, accel_indirectJ2, partials_r, Moon_Z_crs, Sun_Z_crs, GM_Moon,GM_Sun] = force_planets(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config, GM_Earth, radius_Earth, C20)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_planets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbital perturbation due to Sun, Moon and Planets
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


global planets_glob 

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
% Sun, Moon coordinates computation in GCRF & ITRF (m, m/sec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = '3rd_body_peturbations';
[effect_yn] = read_param_cfg(ORB_config,param_keyword);
test_planets = strcmp(effect_yn,'y');

param_keyword = 'Tides_effects';
[effect_yn] = read_param_cfg(ORB_config,param_keyword);
test_tides = strcmp(effect_yn,'y');

if test_planets == 1 || test_tides == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    planets_struct = planets_glob;
    fmod_planets = planets_struct.bodies_list_yn;
    DExxxfilename = planets_struct.ephemeris_filename;
    DEcheby  = planets_struct.DE_chebyshev;
    DErecord = planets_struct.DE_datarecords;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DExxx Reading :
    % 1. DE files (DE data and Header) or
    % 2. global matrices "DEcheby" "DErecord"
    % Select : 1    
    if jd > DErecord(3,1) || jd < DErecord(2,1)
        [DEcheby,DErecord,deformat] = dexxxeph_read(DEfilename,HDfilename,jd);    
        planets_glob.DE_chebyshev = DEcheby;
        planets_glob.DE_datarecords = DErecord;
    end
    % Select : 2
    deln = length(DExxxfilename);
    HDxxxfilename = ['header.' DExxxfilename(deln-2:deln)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sun and Moon Geocentric coordinates
    [rgMoon,vgMoon] = dexxxeph_stategcrf(10,jd,HDxxxfilename,DEcheby,DErecord);
    [rgSun,vgSun] = dexxxeph_stategcrf(11,jd,HDxxxfilename,DEcheby,DErecord);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sun Velocity : Relativistic corrections
    JD2 = jd + 1 / (24*3600);
    [rgSun2,vgSun2] = dexxxeph_stategcrf(11,JD2,HDxxxfilename,DEcheby,DErecord);
    VgSun_dfr = (rgSun2 - rgSun) / 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Body-Fixed positions of Sun and Moon (Transformation : GCRS to ITRS) 
    rgMoon_ITRS = (eopmatrix)' * rgMoon;
    rgSun_ITRS = (eopmatrix)' * rgSun;
    % GM of Sun & Moon in m^3/sec^2 (converted from au^3/d^2 to m^3/sec^2)
    [GMconstant] = dexxxeph_readhd(HDxxxfilename,GM_Earth);
    GMmoon = GMconstant(10,1);
    GMsun  = GMconstant(11,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sun, Moon and Planets 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = '3rd_body_peturbations';
[effect_yn] = read_param_cfg(ORB_config,param_keyword);
test = strcmp(effect_yn,'y');
if test == 1    
    a_bodies = [0; 0; 0];
    pdv_3rdbodies = zeros(3,3);
    for i = 1 : 11
        if fmod_planets(1,i) == 1
            body_ith = i;
            if body_ith == 10
                % Moon
                rg3rd_meter = rgMoon;
            elseif body_ith == 11
                % Sun
                rg3rd_meter = rgSun;
            else                
                % Celestial body's state vector
                [rg3rd_meter,vg3rd_meter] = dexxxeph_stategcrf(body_ith,jd,HDxxxfilename,DEcheby,DErecord);
            end
            % GMbody in m^3/sec^2 (converted from au^3/d^2 to m^3/sec^2)
            [GMconstant] = dexxxeph_readhd(HDxxxfilename);
            GMbody = GMconstant(body_ith,1);
            % Perturbing acceleration
            [abody_x,abody_y,abody_z] = accel_gm3rd(rg3rd_meter,rGCRS,GMbody);
            abody = [abody_x; abody_y; abody_z];
            a_bodies = a_bodies + abody;
            %
            % VEQ : Partial Derivatives
            if VEQ_mode_test == 1
            [pdv3rdbody] = pdv_acclgm3rd(rg3rd_meter,rGCRS,GMbody);
            pdv_3rdbodies = pdv_3rdbodies + pdv3rdbody;
            end
        end
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Indirect J2 effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %C20_J2 = Cnm_GFM(2+1,0+1);
    C20_J2 = C20;
    %[ax_indirectJ2,ay_indirectJ2,az_indirectJ2] = indirectJ2(C20_J2,radius_Earth,GMmoon,rgMoon,GMsun,rgSun);
    [ax_indirectJ2,ay_indirectJ2,az_indirectJ2] = indirectJ2(C20_J2,radius_Earth,GMmoon,rgMoon_ITRS,GMsun,rgSun_ITRS);
    a_indirectJ2 = [ax_indirectJ2 ay_indirectJ2 az_indirectJ2]';
    % Transformation in GCRF 
    a_indirectJ2 = eopmatrix * a_indirectJ2;
    %a_indirectJ2_magn = sqrt(a_indirectJ2(1,1)^2 + a_indirectJ2(2,1)^2 + a_indirectJ2(3,1)^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    a_bodies = [0; 0; 0];
    a_indirectJ2 = [0 0 0]';
    pdv_3rdbodies = zeros(3,3);
end
a_bodies_vec = sqrt(a_bodies(1,1)^2 + a_bodies(2,1)^2 + a_bodies(3,1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accel_bodies = a_bodies;
accel_indirectJ2 = a_indirectJ2;
partials_r = pdv_3rdbodies;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sun and Moon State vectors in Celestial Reference frame
Moon_Z_crs = [rgMoon; vgMoon];
Sun_Z_crs =  [rgSun; vgSun];
% GM constant of Sun and Moon
GM_Moon = GMmoon;
GM_Sun = GMsun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
