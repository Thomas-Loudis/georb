function [orbit_model_struct] = orbit_model_force_parameters (orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model_force_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Forces-related Parameters Number to be estimated per force effect 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                              29 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces Cycle-Per-Revolution (CPR) Periodical Terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
empirical_forces_cpr        = orbit_model_struct.empirical_forces_cpr;
% Number of parameters of empirical forces
Nparam_EMP_FORCE_CPR = empirical_forces_cpr.parameters_number;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer structure array
accelerometer_data_cal_glob = orbit_model_struct.accelerometer_struct;
acc_cal_paramestim_yn = accelerometer_data_cal_glob.param_estim_yn;
test_acc_cal_paramestim = strcmp(acc_cal_paramestim_yn,'y');
if test_acc_cal_paramestim == 1    
% Accelerometer Calibration Parameters matrix
acc_cal_param = accelerometer_data_cal_glob.cal_parameters;
[d1, d2] = size(acc_cal_param);
% Number of Accelerometer calbration model parameters
Nparam_ACC_CAL = d1;
else
Nparam_ACC_CAL = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces (stochstic) Pulses / Piecewise accelerations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulses_stoch_accel_glob     = orbit_model_struct.empirical_forces_pulses;
% Number of pulses / stochastic accelerations' parameters
N_param_pulses_stoch = pulses_stoch_accel_glob.parameters_number;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfm_struct_glob             = orbit_model_struct.gravity_field;
% Number of gravity field parameters (harmonics coefficients) to be estimated
N_param_GRAV = gfm_struct_glob.parameters_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force models/effects with parameters to be estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmodel_PARAM_ESTIM_glob(1) = 0;
Nmodel_PARAM_ESTIM_glob(2) = 0;
Nmodel_PARAM_ESTIM_glob(3) = 0;
Nmodel_PARAM_ESTIM_glob(4) = 0;

% Empirical Forces CPR terms
Nmodel_PARAM_ESTIM_glob(1) = empirical_forces_cpr.effect_01;

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
% Number of parameters to be estimated in addition to the initial state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
