function [obsorbc, obsorbk, obsorbt, obsorbc_ext, obsorbk_ext, obsorbt_ext, obsorbc_full, obsorbk_full, obsorbt_full, COVmatrix] = prm_pseudobs(cfg_fname,GM,eopdat,dpint, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_pseudobs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read/Processing of the orbit data used to form pseudo-observations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_orbit.in :  Orbit modelling configuration file name
% - GM           :  Earth gravity constant  (m^3/sec^2)
%                   Set GM=0 for a standard GM value
% - eopdat       :  Earth Orientation Parameters (EOP) data that are required 
%                   for the orbit arc length
% - dpint        :  Number of data points (days) that are required for the EOP
%                   interpolation to the computation epoch
%
% Output arguments:
% - orbte : State vector at every epoch in ITRS
%           >> z = [MJDt r' v']
% - orbce : State vector at every epoch in GCRS
%           >> z = [MJDt r' v']
% - orbke : Kepler elements at every epoch
%           >> z = [MJDt a e incl Omega omega M]
%
%   position :  r = [x y z]'       in m
%   velocity :  v = [Vx Vy Vz]'    in m/sec
%   epoch    :  MJDt = MJD in TT including the fraction of the day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Set input argument "dpint" equal to 0 for deactivate orbit transformation
% from ITRF to GCRF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                        June  2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 08/07/2022   Thomas Loudis Papanikolaou 
%              Code upgrade 
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2
    % Satellite/Object name    
    param_keyword = 'orbiting_object_name';
    [Object_SAT_ID] = read_param_cfg(cfg_fname,param_keyword);

% Observations modelling
    % Orbit type (kinematic or dynamic)
    param_keyword = 'pseudo_obs_type';
    [pseudo_obs_type] = read_param_cfg(cfg_fname,param_keyword);

    % Pseudo-obs Orbit data file 
    param_keyword = 'pseudo_obs_data';
    [pseudo_obs_data] = read_param_cfg(cfg_fname,param_keyword);
    
    % Covariance matrix of the pseudo-obs orbit data 
    param_keyword = 'cov_pseudo_obs_data';
    [cov_pseudo_obs] = read_param_cfg(cfg_fname,param_keyword);
    
    % Numerical Integrator step    
    param_keyword = 'Stepsize';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_stepsize = sscanf(param_value,'%d %*');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .in configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

    % Satellite/Object name    
    test = strcmp(str1,'orbiting_object_name');
    if test == 1
      Object_SAT_ID = sscanf(line_ith,'%*s %s %*');
    end

% Observations modelling
    % Orbit type (kinematic or dynamic)
    test = strcmp(str1,'pseudo_obs_type');
    if test == 1
      pseudo_obs_type = sscanf(line_ith,'%*s %s %*');
    end

    % Pseudo-obs Orbit data file 
    test = strcmp(str1,'pseudo_obs_data');
    if test == 1
      pseudo_obs_data = sscanf(line_ith,'%*s %s %*');
    end
    
    % Covariance matrix of the pseudo-obs orbit data 
    test = strcmp(str1,'cov_pseudo_obs_data');
    if test == 1
      cov_pseudo_obs = sscanf(line_ith,'%*s %s %*');
    end
    
    % Numerical Integrator step    
    test = strcmp(str1,'Stepsize');
    if test == 1
      integr_stepsize = sscanf(line_ith,'%*s %d %*');
    end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GM == 0
    GM = 3.986004415 * 10^14;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-Observations : Kinematic/Dynamic/Kepler Orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBS data file name
obsorb_fname = pseudo_obs_data; 

% Pseudo-Observations : orbit data type 
% - orbtype : 'dyn' / 'kin' / 'kep'
test = strcmp(pseudo_obs_type,'kinematic');
if test == 1
    orbtype = 'kin';
end
test = strcmp(pseudo_obs_type,'dynamic');
if test == 1
    orbtype = 'dyn';
end
% Kinematic Orbit data format :: Kinematic Orbits by ITSG at TU Graz
test = strcmp(pseudo_obs_type,'itsg_kin');
test_itsg_kin = strcmp(pseudo_obs_type,'itsg_kin');
if test == 1
    %test_itsg_kin = 1;
    orbtype = 'kin';
end
% Kinematic Orbit data format :: GRACE gnv1b
test = strcmp(pseudo_obs_type,'gnv1b');
test_gnv1b = strcmp(pseudo_obs_type,'gnv1b');
if test == 1
    %test_gnv1b = 1;
    orbtype = 'dyn';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obsorb_rate = 1;

% cov_pseudo_obs
test = strcmp(cov_pseudo_obs,'identity');
if test == 1
    COVPform = 0;
end
test = strcmp(cov_pseudo_obs,'diagonal');
if test == 1
    COVPform = 1;
end
test = strcmp(cov_pseudo_obs,'full');
if test == 1
    COVPform = 2;
end

OBS_ORB_Frame = 'trs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematic Orbit read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOCE Kinematic / Reduced-Dynamic orbit
test = strcmp(Object_SAT_ID,'GOCE');
if test == 1    
    tstop = 0;
    [orbte] = sp3c_orb(obsorb_fname,orbtype,tstop);
    if COVPform > 0
        % - COVPform = 1 >> Diagonal
        % - COVPform = 2 >> sub-Full Covariance matrix
        COVdiag = COVPform ;
        SST_PCV2_fname =  ' ' ;
        %[COVmatrix,COVepochs] = goce_SSTPCV2(SST_PCV2_fname,tstop,COVdiag);
    elseif COVPform == 0
        COVmatrix = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO Kinematic or Reduced-Dynamic orbit 
test_grace = 0;    
test = strcmp(Object_SAT_ID,'GRACE-A');
if test == 1    
test_grace = 1;
end
test = strcmp(Object_SAT_ID,'GRACE-B');
if test == 1    
test_grace = 1;    
end

test_gracefo = 0;
test = strcmp(Object_SAT_ID,'GRACE-C');
if test == 1    
test_gracefo = 1;    
end
test = strcmp(Object_SAT_ID,'GRACE-D');
if test == 1    
test_gracefo = 1;    
end

if test_grace == 1 || test_gracefo == 1    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE GNV1B orbits 
    if test_gnv1b == 1 
        [gnv1b,COVmatrix] = grace_gnv1b(obsorb_fname);   
        [sz1 sz2] = size(gnv1b);
        orbte = zeros(sz1,4);
        for i = 1 : sz1
            mjdgps = gnv1b(i,1);
            [t,D,M,Y] = MJD_inv(mjdgps);
            tgps = gnv1b(i,2);
            [tutc,tTT] = time_scales_GPS(tgps,mjdgps);      
            mjdTT = mjdgps + (tTT-tgps)/60/60/24;
            orbte(i,:) = [mjdTT gnv1b(i,3:5)];
            if COVPform > 0
                COVmatrix(i*3-2:i*3,1) = [mjdTT mjdTT mjdTT]';
            end
            clear mjdgps t D M Y tutc tTT mjdTT
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO Kineamtic Orbits    
    elseif test_itsg_kin == 1 
        % Kinematic Orbit data by ITSG at TU-Graz
        [orb_array,orb_xyz] = grace_orb_itsg(obsorb_fname);
        % Orbit reference frame
        OBS_ORB_Frame = 'crs';
        % Time scale transformation: GPS to TT scale
        [sz1, sz2] = size(orb_xyz);
        orbce = zeros(sz1,4);    
        for i = 1 : sz1
            mjdgps = orb_xyz(i,1);
            [t,D,M,Y] = MJD_inv(mjdgps);
            tgps = (mjdgps - floor(mjdgps)) * 24 * 3600;
            % [tutc,tTT] = time_scales_GPS(tgps,mjdgps);
            tTT = tgps + 51.184;
            mjdTT = mjdgps + (tTT-tgps)/60/60/24;
            orbce(i,:) = [mjdTT orb_xyz(i,2:4)];
        end      
        % Covariance matrix is not provided 
        COVmatrix = 0;    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB Orbits    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_georb = strcmp(pseudo_obs_type,'georb');
if test_georb == 1
    % GEORB Orbit data
    [orbit_data_matrix] = read_georb_data(obsorb_fname);
    [sz1, sz2] = size(orbit_data_matrix);
    % Orbit reference frame
    OBS_ORB_Frame = 'crs';
    orbce = zeros(sz1,4);
    % MJD in TT including fraction of the day
    orbce(:,1) = orbit_data_matrix(:,1) + orbit_data_matrix(:,2) / 86400;
    % Orbit X coordinate
    orbce(:,2) = orbit_data_matrix(:,3);
    % Orbit Y coordinate
    orbce(:,3) = orbit_data_matrix(:,4);
    % Orbit Z coordinate
    orbce(:,4) = orbit_data_matrix(:,5);    
    % Covariance matrix is not provided
    COVmatrix = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations Epochs array within the orbit arc length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_fname);
% Orbit model matrix 
IC_MJDo = orbit_model_struct.IC_MJD ;
orbit_arc_length = orbit_model_struct.orbit_arc_length_sec;
% Integration Stepsize
INTGstep = integr_stepsize;
% Orbit arc length in seconds
arc = orbit_arc_length;
% Initial epoch MJD
MJDo = IC_MJDo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBSepochs array
[to_sec,D,M,Y] = MJD_inv(MJDo);
tmax_sec = to_sec + arc;
iobsepoch = 0;
tith_sec = to_sec;
while tith_sec <= tmax_sec
    tith_sec = to_sec + iobsepoch * obsorb_rate * INTGstep;
    iobsepoch = iobsepoch + 1;
end

% Preallocation
OBSepochs = zeros(iobsepoch ,1);

iobsepoch = 0;
tith_sec = to_sec;
while tith_sec <= tmax_sec
    tith_sec = to_sec + iobsepoch * obsorb_rate * INTGstep;
    %OBSepochs(iobsepoch,1) = tith_sec;
    [JDith,MJDith] = MJD_date(tith_sec,D,M,Y);
    OBSepochs(iobsepoch+1,1) = MJDith;
    iobsepoch = iobsepoch + 1;
end
[OBSepochs2] = mjd2mjdtt(OBSepochs,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit transformation from ITRS to GCRS
if OBS_ORB_Frame == 'trs'
    [orbce] = orbt2c(orbte,eopdat,dpint, orbit_model_struct);
    
% Orbit transformation from ICRS to ITRS
elseif OBS_ORB_Frame == 'crs'
    [orbte] = orbc2t(orbce,eopdat,dpint, orbit_model_struct);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d1, d2] = size(orbte);
if d2 < 7
    Vel_vector = 0;
elseif d2 >= 7
    Vel_vector = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Vel_vector == 1
    % Keplerian Elements computations
    [sz1, sz2] = size(orbce);
    % Matrix variable preallocation
    orbke = zeros(sz1,7);
    for ik = 1 : sz1
        [a, e, incl, Omega, omega, f, Mn, E, u] = kepler(orbce(ik,2:4)',orbce(ik,5:7)',GM);
        % Keplerian Elements array
        orbke(ik,:) = [orbce(ik,1) a e incl Omega omega Mn];
    end
elseif Vel_vector == 0
    orbke = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude OBS according to obsorb_rate
obsorbt_full = orbte;
obsorbc_full = orbce;

% Velocity vector
if Vel_vector == 0
    obsorbk_full = obsorbt_full;
else
    obsorbk_full = orbke;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if obsorb_rate == 1
    obsorbc_ext = obsorbc_full;
    obsorbk_ext = obsorbk_full;
    obsorbt_ext = obsorbt_full;
    obsorbc = obsorbc_full;
    obsorbt = obsorbt_full;
    obsorbk = obsorbk_full;
else    
[sz1, sz2] = size(obsorbt_full);
[d1, d2] = size(OBSepochs);
% Preallocation
obsorbt = zeros(d1,sz2);
obsorbc = zeros(d1,sz2);
obsorbk = zeros(d1,sz2);
obsorbt_ext = zeros(sz1-d1,sz2);
obsorbc_ext = zeros(sz1-d1,sz2);
obsorbk_ext = zeros(sz1-d1,sz2);

lp2 = 1;
iobsreduced = 1;
iobs_all0 = 0;
iobs_ext = 0;
iobs = 0;
while iobsreduced <= d1
    MJDobs_reduced = OBSepochs(iobsreduced,1);
    iobs_all = iobs_all0;
    while lp2 == 1 
        iobs_all = iobs_all + 1;
        if iobs_all > sz1
            lp2 = 0;
        else
            MJDobs_all = obsorbt_full(iobs_all,1);
            if abs(MJDobs_all - MJDobs_reduced) < 10^-8
                iobs = iobs + 1;
                obsorbt(iobs,:) = obsorbt_full(iobs_all,:);
                obsorbc(iobs,:) = obsorbc_full(iobs_all,:);
                obsorbk(iobs,:) = obsorbk_full(iobs_all,:);
                lp2 = 0;
                iobs_all0 = iobs_all;
            else
                iobs_ext = iobs_ext + 1;
                obsorbt_ext(iobs_ext,:) = obsorbt_full(iobs_all,:);
                obsorbc_ext(iobs_ext,:) = obsorbc_full(iobs_all,:);
                obsorbk_ext(iobs_ext,:) = obsorbk_full(iobs_all,:);
            end
        end
    end
    iobsreduced = iobsreduced + 1;
    lp2 = 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Vel_vector == 0
    obsorbk_full = 0;
    obsorbk_ext = 0;
    obsorbk = 0;
end
