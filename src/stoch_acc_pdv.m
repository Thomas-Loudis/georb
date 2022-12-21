function [Fpulse, PDr, PDv, PD_param] = stoch_acc_pdv (rsat, vsat, mjd_t, t_sec, delta_v, dir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: pd_pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Acceleration vector and Partial derivatives of pseudo-stochastic pulses 
%  Velocity or Acceleration changes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - rsat:			Satellite Position vector (m)   in inertial frame (ICRF)
% - vsat:			Satellite Velocity vector (m/s) in inertial frame (ICRF)
% - mjd_t:			MJD of current Epoch
% - t_sec:			Seconds since start of day of current Epoch
% - delta_v:		Pulse value
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

global pulses_stoch_accel_glob
%global SCA_array_glob SCA_dpint_glob
global accelerometer_data_cal_glob

pulses_accel_struct = pulses_stoch_accel_glob;

% Matrices intialisation preallocation
PD_param = zeros(3,1);

% State Vector in ICRF
r_crf = rsat;
v_crf = vsat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference frame of pulses axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - 'inertial'      :: Inertial reference frame (X,Y,Z)
% - 'orbital_frame' :: Orbital reference frame (radial,along-track,cross-track)
% - 'srf_frame'     :: Spacecraft reference frame 
PULSES_frame = pulses_accel_struct.reference_frame;

test = strcmp(PULSES_frame,'inertial');
if test == 1
    % Rotational matrix
    %Rmatrix = eye(3,3);
    Rmatrix = [1 0 0; 0 1 0; 0 0 1];    
end

test = strcmp(PULSES_frame,'orbital_frame');
if test == 1
    % Transformation matrix from inertial to orbital frame
    [Rrtn,er,et,en] = orbital_transf(r_crf , v_crf);
    % Transformation matrix from orbital to inertial frame
    Rrtn2crf = inv(Rrtn);
    % Rrtn2crf_T = Rrtn'
    % Rotational matrix
    Rmatrix = Rrtn2crf;   
end

test = strcmp(PULSES_frame,'srf_frame');
if test == 1
    mjd = mjd_t;   
    % Star camera data
    accelerometer_struct = accelerometer_data_cal_glob;
    sca_array = accelerometer_struct.SCA1B_data_array;
    scadpint = accelerometer_struct.sca_interp_no;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Quaternions Interpolation from Star camera data (sca1b)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lagrangian Interpolation
    [qo_int] = interp_Lagrange(sca_array(:,1),sca_array(:,3),mjd,scadpint);
    [q1_int] = interp_Lagrange(sca_array(:,1),sca_array(:,4),mjd,scadpint);
    [q2_int] = interp_Lagrange(sca_array(:,1),sca_array(:,5),mjd,scadpint);
    [q3_int] = interp_Lagrange(sca_array(:,1),sca_array(:,6),mjd,scadpint);
    %sca_int(k,:) = [acc(i,1) acc(i,2) qo_int q1_int q2_int q3_int];
    quaternions_int = [qo_int q1_int q2_int q3_int]';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rotation matrix from SRF (Spacecraft reference frame) to ICRF (Celestial
    % Reference Frame)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rotation matrix based on Quaternion
    %q = sca_interp_matrix(i,3:6)';
    quat = quaternions_int;
    [Rmatrix] = quat_rot(quat);
    %R = R';
    %Rinv1 = R'
    Rinv = inv(Rmatrix);
    Rsrf2crf = Rinv;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
    Rmatrix = Rsrf2crf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit vector of the acceleration axis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dir == 1) 
% Radial direction
delta_v_ti = delta_v(1,1);
e_unit_dir(1,1) = Rmatrix(1,1);
e_unit_dir(2,1) = Rmatrix(2,1);
e_unit_dir(3,1) = Rmatrix(3,1);

elseif (dir == 2) 
% Tangential direction
delta_v_ti = delta_v(2,1);
e_unit_dir(1,1) = Rmatrix(1,2);
e_unit_dir(2,1) = Rmatrix(2,2);
e_unit_dir(3,1) = Rmatrix(3,2);

elseif (dir == 3) 
% Normal direction
delta_v_ti = delta_v(3,1);
e_unit_dir(1,1) = Rmatrix(1,3);
e_unit_dir(2,1) = Rmatrix(2,3);
e_unit_dir(3,1) = Rmatrix(3,3);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirac's "delta" function at the current epoch t (mjd_t)
delta_dirac = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum of Pseudo-stochastic pulses	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fpulse(1,1) = delta_v_ti * delta_dirac * e_unit_dir(1,1);
Fpulse(2,1) = delta_v_ti * delta_dirac * e_unit_dir(2,1);
Fpulse(3,1) = delta_v_ti * delta_dirac * e_unit_dir(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives w.r.t state vector
PDr = 0;
PDv = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives w.r.t unknown parameters to be estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PD_param(1,1) = delta_dirac * e_unit_dir(1,1);
PD_param(2,1) = delta_dirac * e_unit_dir(2,1);
PD_param(3,1) = delta_dirac * e_unit_dir(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
