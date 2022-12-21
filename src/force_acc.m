function [accel_vec, partials_r, partials_p] = force_acc(mjd,Z_crs,Rtrs2crs, EQ_mode, ORB_config)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: force_acc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Accelerometer data processing and calibration modelling for capturing
%  the non-gravitational effects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - mjd:            Epoch's Modified Julian Day number including fraction of the day) in TT (Terrestrial Time)
% - Z_crs:          State vector in Celestial Reference Frame (GCRS)   
%                   z = [r' v']
%                   r:   Position Vector (m) 
%                   v:   Velocity vector (m/sec) 
% - Rtrs2crs:       Tranformation matrix: Terrestrial Reference Frame to Celestial Reference Frame 
% - EQ_mode:        Inegral Equations mode: Equation of Motion or Variational Equations 
% - Configuration array:    Orbit master configuration structure array
%
% Output arguments:
% - accel_vec:      Acceleration vector cartesian components in GCRS
% - partials_r:     Partials w.r.t. position vector 
% - partials_p:     Partials w.r.t. parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Thomas Loudis Papanikolaou                                     
% Created: 15 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Code extracted from function force_eqm_veq. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global accelerometer_data_cal_glob 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations mode: EQM or VEQ
VEQ_mode_test = strcmp(EQ_mode,'VEQ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JD of computation epoch
jd = mjd + 2400000.5;
% Position vector in GCRS
rGCRS = Z_crs(1:3,1);
% Velocity vector in GCRS
vGCRS = Z_crs(4:6,1); 
% Terrestrial (ITRS) to Celestial (GCRS) 
eopmatrix = Rtrs2crs;
% Position vector in ITRS
rITRS = (eopmatrix)' * rGCRS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Satellite Accelerometer data (non-gravitational effects)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
param_keyword = 'acc_data';
[acc_data] = read_param_cfg(ORB_config,param_keyword);
test_accelerometer = strcmp(acc_data,'y');
if test_accelerometer == 1    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHAMP
%     [t,D,M,Y] = MJD_inv(mjd);
%     [facc_x,facc_y,facc_z] = accel_acc(mjd,t,accl,4);
%     facc = [ facc_x; facc_y; facc_z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE and GRACE Follow-On missions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer Calibration Parameters matrix
accelerometer_struct = accelerometer_data_cal_glob;
acc_cal_param = accelerometer_struct.cal_parameters;
acc1b_array_TT = accelerometer_struct.ACC1B_data_array;
acc_dpint = accelerometer_struct.acc_interp_no;
sca1b_array = accelerometer_struct.SCA1B_data_array;
sca_dpint = accelerometer_struct.sca_interp_no;

[PD_ACC_Cal_Param, facc_ICRF, facc_SRF] = pdv_grace_acc(mjd, acc_cal_param, acc1b_array_TT, sca1b_array, acc_dpint, sca_dpint);
pdv_nongrav = zeros(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    facc_ICRF = [0; 0; 0];
    % VEQ pdv
    pdv_nongrav = zeros(3,3); 
    PD_ACC_Cal_Param = zeros(3,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments
accel_vec = facc_ICRF;
partials_r = pdv_nongrav;
partials_p = PD_ACC_Cal_Param;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

