function [SFpulses, PD_pulses_param, PD_pulses_r, PD_pulses_v] = stoch_accel (rsat, vsat, mjd_t, t_sec, stoch_time_interval, pulses_matrix, pulses_axes_vec_01, N_param_pulses, emp_accel_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: stoch_accel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Acceleration vector and Partial derivatives of pseudo-stochastic pulses 
%  Velocity changes or accelerations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - rsat:			Satellite Position vector (m)   in inertial frame (ICRF)
% - vsat:			Satellite Velocity vector (m/s) in inertial frame (ICRF)
% - mjd_t:			MJD of current Epoch
% - t_sec:			Seconds since start of day of current Epoch
% - delta_v:		        Pulse value
% - mjd_ti:			MJD of pulse epoch
% - ti_sec:			Seconds since start of the day of pulse' epoch
% - dir:			Pulse' direction e.g. Radial, Tangential, Normal directions, XYZ directions  
%
% Output arguments:
% - SFpulses:           Acceleration vector cartesian components in inertial frame (GCRF)
% - PD_pulses_param: 	Partial derivatives matrix of the acceleration w.r.t. unknown parameters (GCRF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                           1 September 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Empirical accelerations type: Piecewise constant accelerations or Pulses
% - 'stoch_pulses'          : Pulses as instanteous velocity vector changes at predefined epochs
% - 'stoch_accel_constant'  : Piecewise constant accelerations at predefined epochs 
test = strcmp(emp_accel_type,'stoch_pulses');
if test == 1
emp_accel_id_no = 1;
end
test = strcmp(emp_accel_type,'stoch_accel_constant');
if test == 1
emp_accel_id_no = 2;
end

% Number of Pulses frame' axes conisdered
N_pulses_axes = pulses_axes_vec_01(1) + pulses_axes_vec_01(2) + pulses_axes_vec_01(3);
% Matrices intialisation preallocation
PD_pulses_param = zeros(3 , N_param_pulses);
% Pulses' number of epochs
[N_PULSE_epochs, N2] = size(pulses_matrix);
% Acceleration duration
time_interval = stoch_time_interval;

%PD_pulses_r_sum = zeros(3,3);
%PD_pulses_v_sum = zeros(3,3);
SFpulses(1,1) = 0;
SFpulses(2,1) = 0;
SFpulses(3,1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of loop of empirical accelerations' epochs
for i1 = 1 : N_PULSE_epochs 
   mjd_ti = pulses_matrix(i1,1);
   ti_sec = pulses_matrix(i1,2);

% Conditions control    
    % Dirac's "delta" function at the current epoch t (mjd_t)
    delta_dirac = 0.0D0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulses (instant veclity changes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if emp_accel_id_no == 1
    % Time criteria control 
%     if ( abs(mjd_t - mjd_ti) < 1.0D-06 ) 
%         delta_dirac = 1.0D0;
%     end
    if ( fix(mjd_t) == fix( mjd_ti) ) && (abs(t_sec - ti_sec) < 1.0D-05 )
        delta_dirac = 1.0D0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piece-wise constant accelerations based on step size and time duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif emp_accel_id_no == 2 
    % Time Interval - Acceleration duration  
    dt_accel = time_interval;
    if ( fix(mjd_t) == fix(mjd_ti) ) && (t_sec > ti_sec) && (t_sec <= ti_sec + dt_accel) 
    %if ( fix(mjd_t) == fix(mjd_ti) ) && (t_sec >= ti_sec) && (t_sec < ti_sec + dt_accel) 
    %if ( fix(mjd_t) == fix(mjd_ti) ) && (t_sec >= ti_sec - dt_accel) && (t_sec < ti_sec) 
        delta_dirac = 1.0D0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical acceleration or Pulse vector
if delta_dirac == 1
   % Pulses Vector at epoch ti for the selected directions
   delta_v = zeros(3,1);
   i_dir = 0;
   if ( pulses_axes_vec_01(1) == 1) 
      i_dir = i_dir + 1;
      delta_v(1,1) = pulses_matrix(i1, 2+i_dir);
   end 
   if ( pulses_axes_vec_01(2) == 1)       
      i_dir = i_dir + 1;
      delta_v(2,1) = pulses_matrix(i1 , 2+i_dir);
   end 
   if ( pulses_axes_vec_01(3) == 1) 
      i_dir = i_dir + 1;
      delta_v(3,1) = pulses_matrix(i1 , 2+i_dir);
   end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pulses acceleration per individual axis and overall sum
    i_dir = 0;
    for i2 = 1 : 3
        dir_pulse = i2;        
        if (pulses_axes_vec_01(i2) == 1) 
            i_dir = i_dir + 1;            
            % Call function for computing the stochastic acceleration and its partial derivatives 
            [Fpulse, PD_pulse_r, PD_pulse_v, PD_pulse_param_i] = stoch_acc_pdv (rsat, vsat, mjd_t, t_sec, delta_v, dir_pulse);            
            % Acceleration sum vector 
            SFpulses = SFpulses + Fpulse;
            % Partial derivatives w.r.t. parameters 
            PD_pulses_param(1,i_dir + (i1-1) * N_pulses_axes) = PD_pulse_param_i(1,1);
            PD_pulses_param(2,i_dir + (i1-1) * N_pulses_axes) = PD_pulse_param_i(2,1);
            PD_pulses_param(3,i_dir + (i1-1) * N_pulses_axes) = PD_pulse_param_i(3,1);            
           % Partial derivatives of pulses acceleration w.r.t. state vector
            %PD_pulses_r_sum = PD_pulse_r + PD_pulses_r_sum;
            %PD_pulses_v_sum = PD_pulse_v + PD_pulses_v_sum;
        end
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %break
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Loop of epochs    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PD_pulses_r = PD_pulses_r_sum
%PD_pulses_v = PD_pulses_v_sum
PD_pulses_r = 0;
PD_pulses_v = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
