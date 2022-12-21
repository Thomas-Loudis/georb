function orbit_model (cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit modelling : Read the orbit configuration/parameterisation file and 
%  assign values to global variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 19/04/2021, Thomas Papanikolaou
%             New fucntion global_param.m writen based on code
%             modifications as extracted from the former mainf_DOD.m
%             function (renamed now to orbit_integr.m)  
% 05/07/2022, Thomas Loudis Papanikolaou
%             global_param renamed to orbit_model function. 
%             Major modifications due to the new orbit configuration file format
% 11/12/2022, Thomas Loudis Papanikolaou
%             Code minor modifications reducing global variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelling and Parameterisation data files read | Global variables assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ORB_config_struct_glob
ORB_config_struct_glob = cfg_fname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions and EOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ORBARC_glb MJDo_glb Zo_ICRF_glb
% Read configuration file :: IC and EOP
[orbit_arc_length_sec, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial epoch' Modified Julian Date number
MJDo_glb = IC_MJDo;
% Initial State Vector in ICRF
Zo_ICRF_glb = [IC_MJDo IC_Zo_vec'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit arc length
ORBARC_glb = orbit_arc_length_sec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation Parameters
global EOP_DAT_glob EOP_dpint_glob
EOP_DAT_glob   = EOP_data;
EOP_dpint_glob = EOP_interp_no;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precession-Nutation model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global XYs_IAU200A : defined as global in the function prm_eop.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % List of orbital dynamics
% [gravitational_effects_01, non_gravitational_forces_01, tidal_effects_01, Gravity_field_terms] = prm_dynamics(cfg_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Gravity Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gfm_struct_glob GM_glob
% Gravity model data file (spherical harmonic coefficients)
[GM,ae,Cnm,Snm,sCnm,sSnm, n_max_eqm, m_max_eqm, n_max_veq, m_max_veq, tide_system, gfm_struct] = prm_grav_model(cfg_fname);
gfm_struct_glob = gfm_struct;
GM_glob = GM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planets/Lunar ephemeris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global planets_glob
param_keyword = '3rd_body_peturbations';
[effect_yn] = read_param_cfg(cfg_fname,param_keyword);
test_planets = strcmp(effect_yn,'y');
param_keyword = 'Tides_effects';
[effect_yn] = read_param_cfg(cfg_fname,param_keyword);
test_tides = strcmp(effect_yn,'y');
if test_planets == 1 || test_tides == 1
[planets_struct] = prm_planets(cfg_fname);   
planets_glob = planets_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ocean_tides_struct_glob
param_keyword = 'ocean_tides';
[ocean_tides_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(ocean_tides_yn,'y');
if test == 1
[ocean_tides_struct] = prm_ocean_tides(cfg_fname); 
ocean_tides_struct_glob = ocean_tides_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%--------------------------------------------------------------------------
global aod_struct_glob
param_keyword = 'aod_effects';
[aod_effects_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(aod_effects_yn,'y');
if test == 1 
% Read AOD data
[aod_struct] = prm_aod_data(cfg_fname);
aod_struct_glob = aod_struct;
aod_nmax = aod_struct.degree;
aod_struct_glob.degree_order = [aod_nmax aod_nmax];
aod_struct_glob.degree_eqm = [aod_nmax aod_nmax];
aod_struct_glob.degree_veq = [n_max_veq m_max_veq];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Atmospheric Tides
%--------------------------------------------------------------------------
global aod_tides_struct_glb
param_keyword = 'atm_tides';
[atm_tides_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(atm_tides_yn,'y');
if test == 1
[aod_tides_struct] = prm_aodtides_data(cfg_fname);
aod_tides_struct_glb = aod_tides_struct;
aod_nmax = aod_tides_struct.degree;
aod_tides_struct_glb.degree_order = [aod_nmax aod_nmax];
aod_tides_struct_glb.degree_eqm = [aod_nmax aod_nmax];
aod_tides_struct_glb.degree_veq = [n_max_veq m_max_veq];
end
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global relativistic_effects_glob 
[Relativity_yn, Relativistic_effects_yn, PPN_beta_parameter, PPN_gama_parameter, c_speedoflight] = prm_relativity(cfg_fname);
relativistic_effects_struct.Relativity_yn = Relativity_yn;
relativistic_effects_struct.Relativistic_effects_yn = Relativistic_effects_yn;
relativistic_effects_struct.PPN_beta = PPN_beta_parameter;
relativistic_effects_struct.PPN_gama = PPN_gama_parameter;
relativistic_effects_struct.speedoflight = c_speedoflight;
relativistic_effects_glob = relativistic_effects_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces based on Cycle-Per-Revolution terms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global emp_cpr_glob
[emp_cpr_struct] = prm_empirical_cpr(cfg_fname);
emp_cpr_glob = emp_cpr_struct;
emp_cpr_effect_01 = emp_cpr_struct.effect_01;
% Number of parameters of empirical forces
Nparam_EMP_FORCE_CPR = emp_cpr_struct.parameters_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE mission: Accelerometry Data (Calibration parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global accelerometer_data_cal_glob
param_keyword = 'acc_data';
[acc_data] = read_param_cfg(cfg_fname,param_keyword);
test_accelerometer = strcmp(acc_data,'y');
if test_accelerometer == 1        
% GRACE Accelerometer data preprocessing
[acc_cal_param, acc_dpint, acc1b_array, sca1b_array, sca_dpint, acc_cal_bias_yn, acc_scale_type, accelerometer_struct] = prm_grace_data(cfg_fname);
% Number of accelerometer calibration parameters    
[d1, d2] = size(acc_cal_param);
Nparam_ACC_CAL = d1;
else
    % Accelerometer data use y/n 
    accelerometer_struct.effect_yn = acc_data;
    % Accelerometer calibration parameters estimation y/n 
    accelerometer_struct.param_estim_yn = 'n'; 
    Nparam_ACC_CAL = 0;
end
accelerometer_data_cal_glob = accelerometer_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Accelerations :: Pulses or Piecewise constant accelerations
global pulses_stoch_accel_glob
[pulses_accel_struct] = prm_pulses_stoch(cfg_fname);
pulses_stoch_accel_glob = pulses_accel_struct;
N_param_pulses_stoch = pulses_accel_struct.parameters_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force models/effects with parameters to be estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Nmodel_PARAM_ESTIM_glob
Nmodel_PARAM_ESTIM_glob(1) = 0;
Nmodel_PARAM_ESTIM_glob(2) = 0;
Nmodel_PARAM_ESTIM_glob(3) = 0;

% Empirical Forces CPR terms
Nmodel_PARAM_ESTIM_glob(1) = emp_cpr_effect_01;

% Accelerometer calibration parameters
acc_cal_paramestim_yn = accelerometer_data_cal_glob.param_estim_yn;
test_acc_cal_paramestim = strcmp(acc_cal_paramestim_yn,'y');
if test_acc_cal_paramestim == 1    
    Nmodel_PARAM_ESTIM_glob(2) = 1;
end

% Empirical Accelerations (Piecewise accelerations or Pulses)
PULSES_estim_yn = pulses_stoch_accel_glob.effect_01;
test_empaccel_paramestim = strcmp(PULSES_estim_yn,'y');
if test_empaccel_paramestim == 1    
    Nmodel_PARAM_ESTIM_glob(3) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of unknown parameters to be estimated in addition to the initial
% state vector
global Nparam_GLOB
Nparam_GLOB = 0;

% Empirical Forces
if Nmodel_PARAM_ESTIM_glob(1) == 1
    Nparam_GLOB = Nparam_GLOB + Nparam_EMP_FORCE_CPR;    
end

% Accelerometer calibration parameters
if Nmodel_PARAM_ESTIM_glob(2) == 1
    Nparam_GLOB = Nparam_GLOB + Nparam_ACC_CAL;
end

% Stochastic Pulses parameters
if Nmodel_PARAM_ESTIM_glob(3) == 1
    Nparam_GLOB = Nparam_GLOB + N_param_pulses_stoch;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
