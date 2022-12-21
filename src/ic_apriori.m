function ic_apriori (orbit_config_fname, obsorbc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: ic_apriori.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Initial Conditions: MJD and State Vector obtained from Observations or external orbit data 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in in format 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                 29 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5/11/2022  Dr. Thomas Loudis Papanikolaou
%            Minor changes due to upgrade of the orbit configuration format based on structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
global Zo_ICRF_glb GM_glob   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read configurable parameter :: Initial Conditions
param_keyword = 'ic_state_vector_apriori';
[ic_state_vector_apriori] = read_param_cfg(orbit_config_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Test if Initial Conditions obtained from external orbit 
test_ic_apriori_ext = strcmp(ic_state_vector_apriori,'ext');
% Read External Orbit data based on above conditions 
if test_ic_apriori_ext == 1  
    global EOP_DAT_glob EOP_dpint_glob
    eopdat  = EOP_DAT_glob;
    dpint   = EOP_dpint_glob;
    [orbte,orbce,orbke, ext_orbit_comp_yn] = prm_orbext(orbit_config_fname,GM_glob,eopdat,dpint);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[ic_state_vector_apriori] = read_param_cfg(orbit_config_fname,param_keyword);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions: MJD and State Vector obtained from pseudo-observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if Initial Conditions obtained from Pseudo-Observations (based on kinematic orbit data) 
test_ic_apriori_obs = strcmp(ic_state_vector_apriori,'obs');

% Test if Initial Conditions obtained from external orbit 
test_ic_apriori_ext = strcmp(ic_state_vector_apriori,'ext');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC based on Observations
if test_ic_apriori_obs == 1 || test_ic_apriori_ext == 1
    % Index offset
    i_epoch_to = 1;    
end

if test_ic_apriori_obs == 1
MJD_to_obs = obsorbc(i_epoch_to,1);
r_to_obs   = obsorbc(i_epoch_to,2:4);

% Lagrange Interpolation for computing Velocity vector 
dpint    = 4;
mjd_int  = MJD_to_obs + 1/86400;
[X_int]  = interp_Lagrange(obsorbc(:,1),obsorbc(:,2),mjd_int,dpint);
[Y_int]  = interp_Lagrange(obsorbc(:,1),obsorbc(:,3),mjd_int,dpint);
[Z_int]  = interp_Lagrange(obsorbc(:,1),obsorbc(:,4),mjd_int,dpint);
r_int    = [X_int Y_int Z_int];
v_to_int = r_int - r_to_obs;

% Initial State Vector apriori values
Zo_ICRF_glb(1,1)   = MJD_to_obs;
Zo_ICRF_glb(1,2:4) = r_to_obs;
Zo_ICRF_glb(1,5:7) = v_to_int;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC based on External Orbit 
if test_ic_apriori_ext == 1
% Initial State Vector apriori values based on External Orbit data
Zo_ICRF_glb(1,1)   = orbce(i_epoch_to,1);
Zo_ICRF_glb(1,2:4) = orbce(i_epoch_to,2:4);
Zo_ICRF_glb(1,5:7) = orbce(i_epoch_to,5:7);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC as provided by the input configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_ic_apriori = strcmp(ic_state_vector_apriori,'ic');
if test_ic_apriori == 1
% IC parameters line
param_keyword = 'IC_parameters';
[param_value, param_line] = read_param_cfg(orbit_config_fname,param_keyword);

IC_apriori_rowmatrix = str2num(param_line);
IC_apriori_vecmatrix = IC_apriori_rowmatrix';
ic_apriori_01 = 0;
[Zo_estim, Xaposteriori] = param_aposteriori_apriori(IC_apriori_vecmatrix, ic_apriori_01);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

