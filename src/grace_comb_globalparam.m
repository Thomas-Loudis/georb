function [Nglobal_grav, Nglobal_scale, grav_paramestim_yn] = grace_comb_globalparam(grace_pod_struct, acc_scale_common_01) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function: grace_comb_globalparam 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Number of the global parameters for gravity field and accelerometers scale factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - grace_pod_struct    : GRACE main structure array including orbit model and matrices 
% - acc_scale_common_01 : Accelerometers scale factors flag to consider theses as (1) global/common parameters or (0) arc-based parameters  
% 
% Output arguments: 
% - Nglobal_grav        : Gravity field parameters number (global/common)
% - Nglobal_scale       : Accelerometers scale factors' number of global/common parameters 
% - grav_paramestim_yn  : Gravity parameter estimation flag status y/n on/off 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                          28 September 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%-------------------------------------------------------------------------- 
% Read data from grace_pod structure array 
orbit_model_GRACE1 = grace_pod_struct.grace1_pod.orbit_model; 
orbit_model_GRACE2 = grace_pod_struct.grace2_pod.orbit_model; 
orbit_model_struct = orbit_model_GRACE1;  
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Gravity Field parameters as Global Parameters
%--------------------------------------------------------------------------
% Gravity Field structure array 
gravity_field_struct = orbit_model_struct.gravity_field;
% Gravity solution y/n
gravity_solution_yn = gravity_field_struct.gravity_solution_yn; 
test_status_yn = strcmp(gravity_solution_yn,'y');
COMBESTIM_gravity_param = 0;
if test_status_yn == 1    
COMBESTIM_gravity_param = 1;
end     

N_param_GRAV  = 0;
if COMBESTIM_gravity_param == 1 
grav_paramestim_yn = gravity_field_struct.param_estim_yn; 
test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');    
if test_grav_paramestim_yn == 1    
N_param_GRAV = gravity_field_struct.parameters_number; 
end
else
test_grav_paramestim_yn = 0;
grav_paramestim_yn = 'n'; 
end

% gravity_field_determination_matrix = orbit_model_GRACE1.gravity_field_determination_matrix; 
% COMBESTIM_gravity_step = gravity_field_determination_matrix.COMBESTIM_gravity_step; 
% gfm_struct = orbit_model_GRACE1.gravity_field; 
% N_param_GRAV  = 0;
% if COMBESTIM_gravity_step > 0 
% grav_paramestim_yn = gfm_struct.param_estim_yn; 
% test_grav_paramestim_yn = strcmp(grav_paramestim_yn,'y');    
% if test_grav_paramestim_yn == 1    
% N_param_GRAV = gfm_struct.parameters_number; 
% end
% else
% test_grav_paramestim_yn = 0;
% grav_paramestim_yn = 'n'; 
% end
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Global/Common parameters of Accelerometers Scale factors matrix for both
% GRACE satellites
%-------------------------------------------------------------------------- 
% Accelerometers SCALE matrix :: Commmon/global parameters (overall
% solution e.g. monthly) or non-common/daily parameters (orbit arc solution
% e.g. daily solution)    
% 1. Set to 0 :: non-common (orbit arc solution e.g. daily solution) 
% 2. Set to 1 :: Commmon / Global parameters along with Gravity coefficients solution (e.g. monthly solution)
% N_common_ACC_Scale_01 = 1;
N_common_ACC_Scale_01 = acc_scale_common_01;
%-------------------------------------------------------------------------- 
N_common_ACC_Scale = 0;  
if test_grav_paramestim_yn == 1 && N_common_ACC_Scale_01 == 1  

% Accelerometer Data structure array 
accelerometer_struct = orbit_model_struct.accelerometer_struct; 

% Accelerometer calibration parameters estimation y/n
acc_cal_paramestim_yn = accelerometer_struct.param_estim_yn;  

% Test y/n mode of initial accelerometer calibration parameter estimation
test_acc_cal_paramestim_yn = strcmp(acc_cal_paramestim_yn,'y'); 

if test_acc_cal_paramestim_yn == 1
% Accelerometer Calibration model : Scale factors matrix 
acc_cal_scale_type = accelerometer_struct.cal_scale_type; 

% Scale factors 
scale_matrix_type = acc_cal_scale_type; 
% Diagonal Scale matrix 
SCALE_test = strcmp(scale_matrix_type,'diagonal'); 
if SCALE_test == 1  
    No_acc_scale_param = 3; 
end 
% Semi-Full Scale matrix [Sxy=Syx Sxz=Szx Syz=Szy] 
SCALE_test = strcmp(scale_matrix_type,'semi-full'); 
if SCALE_test == 1  
    No_acc_scale_param = 6; 
end 
% Full Scale matrix 3x3  
SCALE_test = strcmp(scale_matrix_type,'full'); 
if SCALE_test == 1 
    No_acc_scale_param = 9;     
end 

% Number of Scale matrix parameters :: Scale matrix parameters number x2 satellites
% e.g. 3 * 2, 6 * 2, 9 * 2
N_ACC_Scale = 2 * No_acc_scale_param; 

% Number of the Global/Common parameters of Accelerometers Scale factors
N_common_ACC_Scale = N_common_ACC_Scale_01 * N_ACC_Scale;

else
    N_common_ACC_Scale = 0;    
end 

end
%-------------------------------------------------------------------------- 

Nglobal_grav = N_param_GRAV;
Nglobal_scale = N_common_ACC_Scale;



