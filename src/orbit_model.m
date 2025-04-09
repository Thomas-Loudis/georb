function [orbit_model_struct] = orbit_model (cfg_fname,orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit modelling : Read orbit configuration/parameterisation file and
%  forces model data
%  
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
% 04/04/2025, Thomas Loudis Papanikolaou
%             Minor modifications removing all global variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelling and Parameterisation data files read 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration structure array
orbit_model_struct.orbit_config = cfg_fname;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions and EOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read configuration file :: IC and EOP
[orbit_arc_length_sec, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no, IC_Sec_00, TAI_UTC_table, IAU_PN_matrix] = prm_ic(cfg_fname);
% Initial epoch' Modified Julian Date number
orbit_model_struct.IC_MJD = IC_MJDo;
% Initial State Vector in ICRF
orbit_model_struct.IC_CRF = [IC_MJDo IC_Zo_vec'];
% Orbit arc length
orbit_model_struct.orbit_arc_length_sec = orbit_arc_length_sec;
% Earth Orientation Parameters
orbit_model_struct.EOP_data = EOP_data;
orbit_model_struct.EOP_interp_no = EOP_interp_no;
orbit_model_struct.TAI_UTC_table = TAI_UTC_table;
orbit_model_struct.IAU_PN_XYs_matrix = IAU_PN_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precession-Nutation model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XYs_IAU200A : defined within the function prm_eop.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % List of orbital dynamics
% [gravitational_effects_01, non_gravitational_forces_01, tidal_effects_01, Gravity_field_terms] = prm_dynamics(cfg_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Gravity Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity model data file (spherical harmonic coefficients)
[GM,ae,Cnm,Snm,sCnm,sSnm, n_max_eqm, m_max_eqm, n_max_veq, m_max_veq, tide_system, gfm_struct] = prm_grav_model(cfg_fname);
gfm_struct_glob = gfm_struct;
% GM_glob = GM;
N_param_GRAV = gfm_struct.parameters_number;
orbit_model_struct.gravity_field = gfm_struct;
orbit_model_struct.GM_Earth = GM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planets/Lunar ephemeris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = '3rd_body_peturbations';
[effect_yn] = read_param_cfg(cfg_fname,param_keyword);
test_planets = strcmp(effect_yn,'y');
param_keyword = 'Tides_effects';
[effect_yn] = read_param_cfg(cfg_fname,param_keyword);
test_tides = strcmp(effect_yn,'y');
if test_planets == 1 || test_tides == 1
[planets_struct] = prm_planets(cfg_fname);   
% planets_glob = planets_struct;
orbit_model_struct.planets = planets_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_keyword = 'ocean_tides';
[ocean_tides_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(ocean_tides_yn,'y');
if test == 1
[ocean_tides_struct] = prm_ocean_tides(cfg_fname); 
orbit_model_struct.ocean_tides = ocean_tides_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%--------------------------------------------------------------------------
param_keyword = 'aod_effects';
[aod_effects_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(aod_effects_yn,'y');
if test == 1 
% Read AOD data
[aod_struct] = prm_aod_data(cfg_fname);
aod_nmax = aod_struct.degree;
aod_struct.degree_order = [aod_nmax aod_nmax];
aod_struct.degree_eqm = [aod_nmax aod_nmax];
aod_struct.degree_veq = [n_max_veq m_max_veq];
orbit_model_struct.aod_effects = aod_struct;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Atmospheric Tides
%--------------------------------------------------------------------------
param_keyword = 'atm_tides';
[atm_tides_yn] = read_param_cfg(cfg_fname,param_keyword);
test = strcmp(atm_tides_yn,'y');
if test == 1
[aod_tides_struct] = prm_aodtides_data(cfg_fname);
atmTides_nmax = aod_tides_struct.degree;
aod_tides_struct.degree_order = [atmTides_nmax atmTides_nmax];
aod_tides_struct.degree_eqm   = [atmTides_nmax atmTides_nmax];
aod_tides_struct.degree_veq   = [n_max_veq m_max_veq];
orbit_model_struct.atm_tides = aod_tides_struct;
end
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Love numbers (LLN) kn' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LLN_hn, LLN_ln, LLN_kn] = loadlovenumbers();
% Kn' Load Love Numbers
% kn_lln = LLN_kn;
loadlovenumbers_struct.hn = LLN_hn;
loadlovenumbers_struct.ln = LLN_ln;
loadlovenumbers_struct.kn = LLN_kn;
orbit_model_struct.loadlovenumbers = loadlovenumbers_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Pole Tide 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max Degree of spherical harmonics expansion 
nm_values = ocean_tides_struct.degree_eqm;
oceanTides_degree = nm_values(1);
oceanTides_order  = nm_values(2);
opoletide_n = oceanTides_degree;
opoletide_m  = opoletide_n; 
oceanpoletide_struct.opoletide_degree_eqm = [opoletide_n opoletide_m];
oceanpoletide_struct.opoletide_degree_veq = [n_max_veq m_max_veq];

% Ocean Pole Tide :: Coefficients Anm and Bnm including real and imaginary part
% Anm and Bnm coefficients include real and imaginary part' arrays
% AnmR real part
% AnmI Imaginary part
coeffilename = 'desaiscopolecoef.txt';
n_trunc = -1;
[Anm_R,Anm_I,Bnm_R,Bnm_I,Nmax_data] = read_opolecoef(coeffilename,n_trunc); 
oceanpoletide_struct.desaioplecoef_Anm_R = Anm_R;
oceanpoletide_struct.desaioplecoef_Anm_I = Anm_I;
oceanpoletide_struct.desaioplecoef_Bnm_R = Bnm_R;
oceanpoletide_struct.desaioplecoef_Bnm_I = Bnm_I;

orbit_model_struct.oceanpoletide = oceanpoletide_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Relativity_yn, Relativistic_effects_yn, PPN_beta_parameter, PPN_gama_parameter, c_speedoflight] = prm_relativity(cfg_fname);
relativistic_effects_struct.Relativity_yn = Relativity_yn;
relativistic_effects_struct.Relativistic_effects_yn = Relativistic_effects_yn;
relativistic_effects_struct.PPN_beta = PPN_beta_parameter;
relativistic_effects_struct.PPN_gama = PPN_gama_parameter;
relativistic_effects_struct.speedoflight = c_speedoflight;
orbit_model_struct.relativistic_effects = relativistic_effects_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces based on Cycle-Per-Revolution terms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[emp_cpr_struct] = prm_empirical_cpr(cfg_fname);
% emp_cpr_glob = emp_cpr_struct;
emp_cpr_effect_01 = emp_cpr_struct.effect_01;
% Number of parameters of empirical forces
Nparam_EMP_FORCE_CPR = emp_cpr_struct.parameters_number;
orbit_model_struct.empirical_forces_cpr = emp_cpr_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE mission: Accelerometry Data (Calibration parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
orbit_model_struct.accelerometer_struct = accelerometer_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Accelerations :: Pulses or Piecewise constant accelerations
[pulses_accel_struct] = prm_pulses_stoch(cfg_fname);
pulses_stoch_accel_glob = pulses_accel_struct;
N_param_pulses_stoch = pulses_accel_struct.parameters_number;
orbit_model_struct.empirical_forces_pulses = pulses_accel_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force models/effects with parameters to be estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmodel_PARAM_ESTIM_glob(1) = 0;
Nmodel_PARAM_ESTIM_glob(2) = 0;
Nmodel_PARAM_ESTIM_glob(3) = 0;
Nmodel_PARAM_ESTIM_glob(4) = 0;

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

% Gravity Field parameters estimation
GRAV_estim_yn = gfm_struct_glob.param_estim_yn;
test_GRAV_paramestim = strcmp(GRAV_estim_yn,'y');
if test_GRAV_paramestim == 1    
    Nmodel_PARAM_ESTIM_glob(4) = 1;
end

orbit_model_struct.forces_param_estim_yn = Nmodel_PARAM_ESTIM_glob;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of unknown parameters to be estimated in addition to the initial
% state vector
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
% Gravitational parameters
if Nmodel_PARAM_ESTIM_glob(4) == 1
    Nparam_GLOB = Nparam_GLOB + N_param_GRAV;
end
orbit_model_struct.forces_param_estim_no = Nparam_GLOB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
