function [orbte,orbce,orbke, ext_orbit_comp_yn] = prm_orbext(cfg_fname,GM,eopdat,dpint, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read/Processing of the orbit data used for external orbit comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname :  Orbit oncfiguration file (format 2022) 
% - GM        :  Earth gravity constant  (m^3/sec^2)
%                Set GM=0 for a standard GM value
% - eopdat    :  Earth Orientation Parameters (EOP) data that are required 
%                for the orbit arc length
% - dpint     :  Number of data points (days) that are required for the EOP
%                interpolation to the computation epoch
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
% Thomas D. Papanikolaou, AUTH                                  March  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 08/07/2022   Thomas Loudis Papanikolaou 
%              Code upgrade based on the new orbit configuration format
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2    
    % Satellite/Object name    
    param_keyword = 'orbiting_object_name';
    [Object_SAT_ID] = read_param_cfg(cfg_fname,param_keyword);
    
    % Numerical Integrator step    
    param_keyword = 'Stepsize';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_stepsize = sscanf(param_value,'%d %*');    
    
% External orbit comparison
    % Option for applying external orbit comparison
    param_keyword = 'external_orbit_comp';
    [ext_orbit_comp_yn] = read_param_cfg(cfg_fname,param_keyword);

    % External orbit type  ('dynamic' or 'kinematic')  
    param_keyword = 'external_orbit_type';
    [ext_orbit_type] = read_param_cfg(cfg_fname,param_keyword);

    % External orbit data file name     
    param_keyword = 'external_orbit_data';
    [ext_orbit_data_fname] = read_param_cfg(cfg_fname,param_keyword);
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
    
    % Numerical Integrator step    
    test = strcmp(str1,'Stepsize');
    if test == 1
      integr_stepsize = sscanf(line_ith,'%*s %d %*');
    end
    
% External orbit comparison
    % Option for applying external orbit comparison
    test = strcmp(str1,'external_orbit_comp');
    if test == 1
      ext_orbit_comp_yn = sscanf(line_ith,'%*s %s %*');
    end

    % External orbit type  ('dynamic' or 'kinematic')  
    test = strcmp(str1,'external_orbit_type');
    if test == 1
      ext_orbit_type = sscanf(line_ith,'%*s %s %*');
    end

    % External orbit data file name     
    test = strcmp(str1,'external_orbit_data');
    if test == 1
      ext_orbit_data_fname = sscanf(line_ith,'%*s %s %*');
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
% Orbiting Object/Satellite name to ID number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE satellite mission
test_grace = 0;    
test = strcmp(Object_SAT_ID,'GRACE-A');
if test == 1    
test_grace = 1;
end
test = strcmp(Object_SAT_ID,'GRACE-B');
if test == 1    
test_grace = 1;    
end

% GRACE-FO satellite mission
test_gracefo = 0;
test = strcmp(Object_SAT_ID,'GRACE-C');
if test == 1    
test_gracefo = 1;    
end
test = strcmp(Object_SAT_ID,'GRACE-D');
if test == 1    
test_gracefo = 1;    
end

% GOCE satellite mission
test = strcmp(Object_SAT_ID,'GOCE');
if test == 1    
test_groce = 1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit type
test = strcmp(ext_orbit_type,'kinematic');
if test == 1
    ext_orbtype = 'kin';
end

test = strcmp(ext_orbit_type,'dynamic');
if test == 1
    ext_orbtype = 'dyn';
end

% External Orbit data format :: GRACE gnv1b
test_gnv1b = strcmp(ext_orbit_type,'gnv1b');
test_kepler = strcmp(ext_orbit_type,'Kepler'); 
test_goce = strcmp(ext_orbit_type,'goce_pso'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit data read :: Orbit arrays: orbke_o, orbke, orbce, orbte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal orbit (Keplerian Orbit)
if test_kepler == 1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Orbit Ephemeris read
elseif test_goce == 1
    % GOCE Reduced-Dynamic orbit
    orbfname = ext_orbit_data_fname;
    orbtype = ext_orbtype;
    % arc = orbit_arc_length;
    tstop = 0;
    [orbte] = sp3c_orb(orbfname,orbtype,tstop);
    
elseif test_gnv1b == 1
    % GRACE-FO and GRACE GNV1B orbts        
    gnv1bfname = ext_orbit_data_fname;
    [gnv1b] = grace_gnv1b(gnv1bfname);
    [sz1 sz2] = size(gnv1b);
    orbte = zeros(sz1,sz2-1);
    for i = 1 : sz1
        mjdgps = gnv1b(i,1);
        [t,D,M,Y] = MJD_inv(mjdgps);
        tgps = gnv1b(i,2);
        % [tutc,tTT] = time_scales_GPS(tgps,mjdgps);
        tTT = tgps + 51.184;
        mjdTT = mjdgps + (tTT-tgps)/60/60/24;
        orbte(i,:) = [mjdTT gnv1b(i,3:end)];
    end    
end

% Orbit data format
test = strcmp(ext_orbit_type,'itsg_redyn');
if test == 1
    % test_itsg_redyn = 1;
end

test = strcmp(ext_orbit_type,'georb');
if test == 1
    % GEORB orbit data format
    [orbit_data_matrix] = read_georb_data(ext_orbit_data_fname);
    orbte(:,1)   = orbit_data_matrix(:,1) + orbit_data_matrix(:,2) / 86400;
    orbte(:,2:7) = orbit_data_matrix(:,3:8);    
    % % Orbit transformation from GCRS to ITRS
    % orbce = orbte;
    % [orbte] = orbc2t(orbce,eopdat,dpint);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dpint > 0  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit transformation from ITRS to GCRS
[orbce] = orbt2c(orbte,eopdat,dpint,orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keplerian Elements computations
[sz1 sz2] = size(orbce);
for ik = 1 : sz1
    [a(ik,1),e(ik,1),incl(ik,1),Omega(ik,1),omega(ik,1),f(ik,1),Mn(ik,1),E(ik,1),u(ik,1)] = kepler(orbce(ik,2:4)',orbce(ik,5:7)',GM);
end
% Keplerian Elements array
orbke = [orbce(:,1) a e incl Omega omega Mn];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    orbce = 0;
    orbke = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
