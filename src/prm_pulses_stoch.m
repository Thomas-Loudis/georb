function [pulses_accel_struct] = prm_pulses_stoch (cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: cfg_empirical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read configuration file for Empirical Forces variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                               1 September 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
% 14/12/2022  Dr. Thomas Loudis Papanikolaou
%             Code upgrade to pass output arguments to structure array 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2
% Stochastic Pulses 
    param_keyword = 'PULSES_estim_yn';
    [PULSES_estim_yn] = read_param_cfg(cfg_fname,param_keyword);

% Stochastic parameters type (Pulses or piecewise accelerations)    
    param_keyword = 'stoch_param_type';
    [stoch_param_type] = read_param_cfg(cfg_fname,param_keyword);
    
% Pulses orbital frame    
    param_keyword = 'PULSES_frame';
    [PULSES_frame] = read_param_cfg(cfg_fname,param_keyword);

% Number of pulses epochs 
    param_keyword = 'PULSES_epochs_number';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_epochs_number = sscanf(param_value,'%d %*');
        
% PULSES epochs interval 
    param_keyword = 'PULSES_interval';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_interval = sscanf(param_value,'%d %*');

% Time interval of the acceleration (duration) 
    param_keyword = 'stoch_time_interval';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    stoch_time_interval = sscanf(param_value,'%d %*');
    
% Time offset of first and last pulses epochs (seconds)
    param_keyword = 'PULSES_offset';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_offset = sscanf(param_value,'%d %*');
    
% Axes components included (e.g. radial/tangential/normal, X/Y/Z)
% Axis 1 (radial, X)    
    param_keyword = 'PULSES_axis_1';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_axis_1 = sscanf(param_value,'%d %*');
% Axis 2 (radial, X)    
    param_keyword = 'PULSES_axis_2';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_axis_2 = sscanf(param_value,'%d %*');
% Axis 3 (radial, X)    
    param_keyword = 'PULSES_axis_3';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PULSES_axis_3 = sscanf(param_value,'%d %*');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

% Stochastic Pulses 
    test = strcmp(str1,'PULSES_estim_yn');
    if test == 1
      PULSES_estim_yn = sscanf(line_ith,'%*s %s %*') ;
    end

% Stochastic parameters type (Pulses or piecewise accelerations)    
    test = strcmp(str1,'stoch_param_type');
    if test == 1
      stoch_param_type = sscanf(line_ith,'%*s %s %*') ;
    end
    
% Pulses orbital frame    
    test = strcmp(str1,'PULSES_frame');
    if test == 1
      PULSES_frame = sscanf(line_ith,'%*s %s %*') ;
    end

% Number of pulses epochs 
    test = strcmp(str1,'PULSES_epochs_number');
    if test == 1
      PULSES_epochs_number = sscanf(line_ith,'%*s %d %*') ;
    end
    
% PULSES epochs interval 
    test = strcmp(str1,'PULSES_interval');
    if test == 1
      PULSES_interval = sscanf(line_ith,'%*s %d %*') ;
    end

% Time interval of the acceleration (duration) 
    test = strcmp(str1,'stoch_time_interval');
    if test == 1
      stoch_time_interval = sscanf(line_ith,'%*s %d %*') ;
    end
    
% Time offset of first and last pulses epochs (seconds)
    test = strcmp(str1,'PULSES_offset');
    if test == 1
      PULSES_offset = sscanf(line_ith,'%*s %d %*') ;
    end
    
% Axes components included (e.g. radial/tangential/normal, X/Y/Z)
% Axis 1 (radial, X)    
    test = strcmp(str1,'PULSES_axis_1');
    if test == 1
      PULSES_axis_1 = sscanf(line_ith,'%*s %d %*') ;
    end
% Axis 2 (radial, X)    
    test = strcmp(str1,'PULSES_axis_2');
    if test == 1
      PULSES_axis_2 = sscanf(line_ith,'%*s %d %*') ;
    end
% Axis 3 (radial, X)    
    test = strcmp(str1,'PULSES_axis_3');
    if test == 1
      PULSES_axis_3 = sscanf(line_ith,'%*s %d %*') ;
    end
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulses axes vector
pulses_axes_vec_01(1,1) = PULSES_axis_1;
pulses_axes_vec_01(2,1) = PULSES_axis_2;
pulses_axes_vec_01(3,1) = PULSES_axis_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC Initial Epoch
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no, IC_Sec_00] = prm_ic(cfg_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N_pulses_epochs = PULSES_epochs_number;
pulse_step      = PULSES_interval;
pulse_offset    = PULSES_offset;
N_pulses_epochs = fix(orbit_arc_length/pulse_step)+1 - (ceil(pulse_offset/pulse_step));
if N_pulses_epochs <= 0 
    N_pulses_epochs = 1;
end
mjd_epoch0 = IC_MJDo;
sec_epoch0 = IC_Sec_00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if stochastic pulses are included in the orbit modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True
test = strcmp(PULSES_estim_yn,'y');
if test == 1
% Form the Pulses matrix    
    [pulses_matrix_init, N_pulses_axes, N_pulses_param] = pulses_init(N_pulses_epochs, pulses_axes_vec_01, mjd_epoch0, sec_epoch0, pulse_step, pulse_offset);
end

% False
test = strcmp(PULSES_estim_yn,'n');
if test == 1
    N_pulses_param     = 0;
    pulses_matrix_init = 0;
    N_pulses_axes = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force - Effect: y/n 
pulses_accel_struct.effect_01 = PULSES_estim_yn;

% Reference Frame of the empirical acclerations
pulses_accel_struct.reference_frame = PULSES_frame;

% Array for defining the parameters to be estimated
pulses_accel_struct.parameters_01 = pulses_axes_vec_01;

% Empirical accelerations matrix 
pulses_accel_struct.acceleration_matrix = pulses_matrix_init;

% Number of parameters of empirical accelerations/pulses
pulses_accel_struct.parameters_number = N_pulses_param;

% Empirical accelerations type: Piecewise constant accelerations or Pulses
pulses_accel_struct.parameters_type = stoch_param_type;

% Empirical accelerations: Time duration interval 
pulses_accel_struct.step = pulse_step;

% Empirical accelerations: Time duration interval 
pulses_accel_struct.duration = stoch_time_interval;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
