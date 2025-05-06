function [config_struct,orbit_model_struct] = config_struct_update(config_struct, orbit_pod_mode, orbit_parameters, ic_apriori_01, acc_cal_update, grav_param_01, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  config_struct_update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Modify configuration parameters of main configuration structure  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - config_struct     : Configuration structure name 
% - ic_apriori_01 :
%   0 : Apriori values of IC Parameters are provided by the variable Xmatrix
%   1 : Corrections to IC Parameters are provided by the variable Xmatrix
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                              14 April 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
% 29/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit/Force model matrix
accelerometer_data_cal_glob = orbit_model_struct.accelerometer_struct; 
IC_MJDo = orbit_model_struct.IC_MJD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify Variables of the orbit_model and orbit_configuration structure arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity field parameters estimation :: Change mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if grav_param_01 == 1 
grav_field_paramestim_yn = 'y';

% Update Configuration structure array
% Modify :: PARAM grav_field_paramestim_yn :: VALUE y
param_keyword = 'grav_field_paramestim_yn';
param_value = grav_field_paramestim_yn;
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Update orbit model structure array
[orbit_model_struct] = orbit_model_gravity (grav_field_paramestim_yn,orbit_model_struct);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set ACC calibration parameters to fixed values; Cancelled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if acc_cal_update == 100
% Modify :: PARAM acc_cal_paramestim :: VALUE n
param_keyword = 'acc_cal_paramestim';
param_value = 'n';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% ic_apriori_01 = 0;
[Zo_estim, Xaposteriori,orbit_model_struct] = param_aposteriori_apriori(orbit_parameters, ic_apriori_01,orbit_model_struct);

accelerometer_struct = accelerometer_data_cal_glob;

% Accelerometer Calibration Parameters matrix
acc_cal_parameters_ic = accelerometer_struct.cal_parameters;

% param_value :: acc_cal_parameters_ic
param_value   = sprintf('%-27.17e ', acc_cal_parameters_ic);
% Modify :: PARAM keyword
param_keyword = 'acc_cal_parameters_ic';

[k m] = size(config_struct);
i_struct = m;
i_struct = i_struct + 1;
param_keyword = 'acc_cal_parameters_ic';
config_struct(i_struct).names  = param_keyword;
config_struct(i_struct).values = param_value;

% [config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orbit_pod_mode = 'orbit_propagation_veq';
orbit_pod_mode_input = orbit_pod_mode;

% Modify :: PARAM orbit_mode :: VALUE orbit_propagation_veq
param_keyword = 'orbit_mode';
param_value = orbit_pod_mode_input;
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Modify :: PARAM ic_state_vector_apriori :: VALUE ic :: Set to 'ic' responds to be read by the input IC vector
param_keyword = 'ic_state_vector_apriori';
param_value = 'ic';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Modify :: PARAM Reference_frame :: VALUE ICRF
param_keyword = 'Reference_frame';
param_value   = 'ICRF';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Modify :: PARAM Time_scale :: VALUE TT
param_keyword = 'Time_scale';
param_value   = 'TT';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Modify :: PARAM Date_format :: VALUE MJD
param_keyword = 'Date_format';
param_value   = 'MJD';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orbit.model orbit.config update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic_apriori_01 > -1
    
Xaposteriori = orbit_parameters;
if length(Xaposteriori) > 1 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update :: Orbit model matrix 
% ic_apriori_01 = 0;
[Zo_estim, Xaposteriori_2, orbit_model_struct] = param_aposteriori_apriori(orbit_parameters, ic_apriori_01, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ic_apriori_01 == 1
% Xaposteriori = orbit_parameters + Xapriori;
Xaposteriori = Xaposteriori_2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update :: Orbit configuration matrix

MJD_0_fraction = IC_MJDo;  % orbit_matrix(1,1); 
MJD_0  = fix(MJD_0_fraction);
Sec_00 = (MJD_0_fraction - fix(MJD_0_fraction)) * 86400;

% Modify :: PARAM Initial_Epoch :: VALUE 
param_keyword = 'Initial_Epoch';
%param_value   = sprintf('%s ','0','0',Sec_00,MJD_0);
param_value   = sprintf('%d %d %.17f %d',0,0,Sec_00,MJD_0);
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% param_value :: Xaposteriori
param_value   = sprintf('%-27.17e ', Xaposteriori);

% Modify :: PARAM State_vector :: VALUE Xaposteriori
param_keyword = 'State_vector';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);

% Modify :: PARAM IC_parameters :: VALUE Xaposteriori
param_keyword = 'IC_parameters';
[config_struct] = write_configstruct_cor(config_struct, param_keyword, param_value);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orbit_model_struct.orbit_config = config_struct; 
