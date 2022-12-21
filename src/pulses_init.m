function [pulses_matrix_init, N_pulses_axes, N_pulses_param] = pulses_init(N_pulses_epochs, pulses_axes_vec_01, mjd_epoch0, sec_epoch0, pulse_step, pulse_offset)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: pulses_init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Piecewise accelerations or Velocity changes apriori values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - main_config_fname     : GEORB master configuration file name
% - ic_data_object        : 
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                           1 September 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulses axes considered
N_pulses_axes = pulses_axes_vec_01(1,1) + pulses_axes_vec_01(2,1) + pulses_axes_vec_01(3,1);
% Number of pulses to be estimated
N_pulses_param = N_pulses_axes * N_pulses_epochs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrices intialisation preallocation
pulses_matrix_init = zeros(N_pulses_epochs, 2 + N_pulses_axes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulses matrix initialisation to apriori values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1_pulse = 1 : N_pulses_epochs
    if (i1_pulse == 1)        
        offset = pulse_offset;
    elseif (i1_pulse == N_pulses_epochs)
       %offset = -1.0D0 * pulse_offset;
       %offset = pulse_offset; 
       %offset = 0.0D0;
       piecewise_time = 7*60;
       offset = - 15*60; %-2 * piecewise_time;
    else
       offset = 0.0D0;
       %offset = pulse_offset;
    end
    offset = pulse_offset;    
    pulses_matrix_init(i1_pulse, 1) = mjd_epoch0 + offset / 86400 + (i1_pulse-1) * (pulse_step / 86400);
    pulses_matrix_init(i1_pulse, 2) = sec_epoch0 + offset + (i1_pulse-1) * pulse_step;
    for i2_pulse = 1 : N_pulses_axes
        pulses_matrix_init(i1_pulse,i2_pulse+2) = 1.0D-12;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulses_matrix = pulses_matrix_init;
[n1 n2] = size(pulses_matrix);
pulses_parameters_matrix = zeros(N_pulses_param,1);
k = 0;
for i = 1 : n1
    for j = 3 : n2
        k = k + 1;
        pulses_parameters_matrix(k,1) = pulses_matrix(i,j);
    end
end
% Apriori parametes values
pulses_param_apriori = pulses_parameters_matrix;
