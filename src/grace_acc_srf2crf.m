function [acc_array_icrf,acc_array_srf] = grace_acc_srf2crf(acc1b_array, acc_cal_param, sca1b_interp_array)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_acc_preproc:  GRACE accelerometry data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Processing of GRACE accelerometry data by applying calibration
%  parameters and rotation to inertial frame using star camera data. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - grace:    GRACE satellite ID
%             GRACE-A: grace=1, GRACE-B: grace=2 
% - acc:      accelerometry data in SRF (Science Reference Frame)
% - sca:      star camera assembly data (quatenions or Euler angles)
% - scadpint: Number of data points required for "sca" interpolation
% - cal:      calibration parameters (ISDC or DEOS parameters)
%             ISDC : cal = 1
%             DEOS : cal = DEOS_calprmdata_array
% - caldpint: Number of data points required for "cal" interpolation
%
% Output arguments:
% - facc_array:   Calibrated accelerometry array in inertial ICRS
%   facc_array = [MJD_TT tTT lin_accl_x lin_accl_y lin_accl_z]
%   MJDgps:      MJD in TT time scale including fraction of the day
%   tTT:         Seconds since 0h in TT time scale
%   lin_accl_i:  Calibrated linear acceleration along ICRS axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 22/04/2021, Dr. Thomas Papanikolaou
%             Renamed and upgraded from former grace_accproc.m function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial Accelerometry data array
% acc rate 1 sec
[sz1 sz2] = size(acc1b_array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation from SRF to ICRS and implementation of Calibration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc_array_icrf    = zeros(sz1,5);
acc_array_srf = zeros(sz1,5);
for i = 1 : sz1
% Convert time scale from GPS time to Terrestrial Time scale
    [tgps,D,M,Y] = MJD_inv(acc1b_array(i,1));
    tgps = acc1b_array(i,2);
    [tutc,tTT] = time_scales_GPS(tgps,acc1b_array(i,1));
    mjdTT = acc1b_array(i,1) + (tTT - tgps) / (24*3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Calibration parameters (biases & scales) to accelerometry array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    facc_srf = [mjdTT tTT (acc_cal_param(1,1)+acc_cal_param(4,1) *acc1b_array(i,3)) ... 
                          (acc_cal_param(2,1)+acc_cal_param(5,1) *acc1b_array(i,4)) ... 
                          (acc_cal_param(3,1)+acc_cal_param(6,1) *acc1b_array(i,5))];
    acc_array_srf(i,:) = facc_srf;                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation from SRF to ICRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply rotation parameters (star camera data) to the calibrated
% accelerometry array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rotation matrix based on Quaternion
    q = sca1b_interp_array(i,3:6)';
    [R] = quat_rot(q);
    Rinv = inv(R);
    
    f_srf = facc_srf(1,3:5)';
    f_icrs = Rinv * f_srf;
    mjdTT = facc_srf(1,1);
    tTT = facc_srf(1,2);
    acc_array_icrf(i,:) = [mjdTT tTT f_icrs(1,1) f_icrs(2,1) f_icrs(3,1)];
end
