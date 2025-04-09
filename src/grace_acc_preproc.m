function [acc_cal_param, acc1b_array_TT, sca1b_array_TT] = grace_acc_preproc(acc1b_array,sca1b_array,acc_cal_model,scale_matrix_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_acc_preproc:  GRACE accelerometry data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Processing of GRACE accelerometry data by applying calibration
%  parameters and rotation to inertial frame using star camera data. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - acc:                accelerometer data in SRF (Science Reference Frame)
% - sca:                star camera assembly data (quatenions or Euler angles)
% - acc_cal_model:      
% - scale_matrix_type:   
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
% GRACE Accelerometry Calibration Parameters: Initial values (apriori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[acc_cal_param] = grace_acc_cal_init(acc1b_array,acc_cal_model,scale_matrix_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time system transformation from GPS Time to Terrestrial Time (TT) of the
% initial arrays of acc1b (accelerometry data) and sca1b (star camera data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acc1b array
[sz1 sz2] = size(acc1b_array);
acc1b_array_TT = zeros(sz1,sz2);
for i = 1 : sz1
% Convert time scale from GPS time to Terrestrial Time scale
    tgps = acc1b_array(i,2);
    % Difference between Terrestrial Time and GPS time scale in sec
    dt_TT_GPS_sec = 51.184;
    mjdTT = acc1b_array(i,1) + dt_TT_GPS_sec / (24*3600);
    tTT = tgps + dt_TT_GPS_sec;
    %[JD,MJD] = MJD_date(t,D,M,Y)
    acc1b_array_TT(i,:) = [mjdTT tTT acc1b_array(i,3) acc1b_array(i,4) acc1b_array(i,5)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sca1b_array 
[sz1 sz2] = size(sca1b_array);
sca1b_array_TT = zeros(sz1,sz2);
for i = 1 : sz1
% Convert time scale from GPS time to Terrestrial Time scale
    tgps = sca1b_array(i,2);
    dt_TT_GPS_sec = 51.184;
    mjdTT = acc1b_array(i,1) + dt_TT_GPS_sec / (24*3600);
    tTT = tgps + dt_TT_GPS_sec;
    %[JD,MJD] = MJD_date(t,D,M,Y)
    sca1b_array_TT(i,:) = [mjdTT tTT sca1b_array(i,3) sca1b_array(i,4) sca1b_array(i,5) sca1b_array(i,6)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
