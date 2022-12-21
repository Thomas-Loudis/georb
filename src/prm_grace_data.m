function [acc_cal_param, acc_dpint, acc1b_array_TT, sca1b_array_TT, sca_dpint, acc_cal_bias_yn, acc_cal_scale_type, accelerometer_struct] = prm_grace_data(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  prm_grace_data.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
% Reading and setting of the parameters/variables required for the GRACE
% data preprocessing, in particualr the GRACE accelerometry data
% calibration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename : prm.in input file name
%
% Output arguments:
% - Fmod    : Dynamic effects
%             >> Fmod = [prmfmod(1,1); prmfmod(1,2); prmfmod(1,3);
%                        prmfmod(1,4); prmfmod(1,5); prmfmod(1,6)]; 
% - FmodEGM : EGM parameters
%             >> FmodEGM = [n m GM ae egmid]';
% - Cnm     : Spherical Harmonic Coefficients (SHC) array
% - Snm     : Spherical Harmonic Coefficients (SHC) array
% - sCnm    : Variances array of Cnm
% - sSnm    : Variances array of Snm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Papanikolaou                                     22 April 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 22/04/2021, Dr. Thomas Papanikolaou
%             Code extracted from prm3.m function into an individual function
% 06/07/2022, Dr. Thomas Loudis Papanikolaou
%             Code update to read the new format of the orbit configuration
%             file
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2
% Satellite ID name 
    param_keyword = 'orbiting_object_name';
    [grace_sat_id] = read_param_cfg(cfg_fname,param_keyword);
    
% Accelerometer data use y/n 
    param_keyword = 'acc_data';
    [acc_data_yn] = read_param_cfg(cfg_fname,param_keyword);
    
% Accelerometer calibration parameters estimation y/n 
    param_keyword = 'acc_cal_paramestim';
    [acc_cal_paramestim_yn] = read_param_cfg(cfg_fname,param_keyword);
    
% Accelerometer data file name
    param_keyword = 'accelerometer_data';
    [acc_data_fname] = read_param_cfg(cfg_fname,param_keyword);

% Accelerometry data interpolator number of data points 
    param_keyword = 'accelerometer_interp_no';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    acc_dpint = sscanf(param_value,'%d %*');
    
% Star Camera data 
    param_keyword = 'star_camera_data';
    [sca_data_fname] = read_param_cfg(cfg_fname,param_keyword);

% Star Camera data interpolator number of data points 
    param_keyword = 'star_camera_interp_no';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    sca_dpint = sscanf(param_value,'%d %*');
    
% Accelerometer Calibration parameters 
%  Bias parameters
    param_keyword = 'acc_cal_bias';
    [acc_cal_bias] = read_param_cfg(cfg_fname,param_keyword);
    
% Bias drift 1st order
    param_keyword = 'acc_cal_bias_drift_1';
    [acc_cal_bias_drift_1] = read_param_cfg(cfg_fname,param_keyword);
    
% Bias drift 2nd order
    param_keyword = 'acc_cal_bias_drift_2';
    [acc_cal_bias_drift_2] = read_param_cfg(cfg_fname,param_keyword);

% Scale matrix
    param_keyword = 'acc_cal_scale';
    [acc_cal_scale] = read_param_cfg(cfg_fname,param_keyword);
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

% Satellite ID name 
    test = strcmp(str1,'orbiting_object_name');
    if test == 1
      grace_sat_id = sscanf(line_ith,'%*s %s %*') ;
    end
    
% Accelerometer data 
    test = strcmp(str1,'accelerometer_data');
    if test == 1
      acc_data_fname = sscanf(line_ith,'%*s %s %*') ;
    end

% Accelerometry data interpolator number of data points 
    test = strcmp(str1,'accelerometer_interp_no');
    if test == 1
      acc_dpint = sscanf(line_ith,'%*s %d %*');
    end
    
% Star Camera data 
    test = strcmp(str1,'star_camera_data');
    if test == 1
      sca_data_fname = sscanf(line_ith,'%*s %s %*') ;
    end

% Star Camera data interpolator number of data points 
    test = strcmp(str1,'star_camera_interp_no');
    if test == 1
      sca_dpint = sscanf(line_ith,'%*s %d %*');
    end
    
    
% Accelerometer Calibration parameters 
%  Bias parameters
    test = strcmp(str1,'acc_cal_bias');
    if test == 1
      acc_cal_bias = sscanf(line_ith,'%*s %s %*') ;
    end
    
% Bias drift 1st order
    test = strcmp(str1,'acc_cal_bias_drift_1');
    if test == 1
      acc_cal_bias_drift_1 = sscanf(line_ith,'%*s %s %*') ;
    end
    
% Bias drift 2nd order
    test = strcmp(str1,'acc_cal_bias_drift_2');
    if test == 1
      acc_cal_bias_drift_2 = sscanf(line_ith,'%*s %s %*') ;
    end

% Sace matrix
    test = strcmp(str1,'acc_cal_scale');
    if test == 1
      acc_cal_scale = sscanf(line_ith,'%*s %s %*') ;
    end   
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleremetry Calibration modelling auxiliary matrices : Bias parameters
test = strcmp(acc_cal_bias,'y');
acc_cal_bias_yn(1) = test;

test = strcmp(acc_cal_bias_drift_1,'y');
acc_cal_bias_yn(2) = test;

test = strcmp(acc_cal_bias_drift_2,'y');
acc_cal_bias_yn(3) = test;

acc_cal_scale_type = acc_cal_scale;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometry (Satellite on-board accelerometer measurements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstop = 0;

% Read Accelerometry 1b data
[acc1b_array] = grace_acc1b(acc_data_fname,tstop);

% Read Star Camera 1b data
[sca1b_array] = grace_sca1b(sca_data_fname,tstop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometry data preprocessing and calibration parameters initialisation
[acc_cal_param, acc1b_array_TT, sca1b_array_TT] = grace_acc_preproc(acc1b_array,sca1b_array,acc_cal_bias_yn,acc_cal_scale_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer data use y/n 
accelerometer_struct.effect_yn = acc_data_yn;
% Accelerometer calibration parameters estimation y/n 
accelerometer_struct.param_estim_yn = acc_cal_paramestim_yn;
% Accelerometer Calibration Parameters matrix
accelerometer_struct.cal_parameters = acc_cal_param;
% Accelerometer Calibration model parameters to be included and estimated
accelerometer_struct.cal_bias_yn = acc_cal_bias_yn;
accelerometer_struct.cal_scale_type = acc_cal_scale_type;
% Data arrays
accelerometer_struct.ACC1B_data_array = acc1b_array_TT;
accelerometer_struct.SCA1B_data_array = sca1b_array_TT;
% Lagrange Interpolation polynomial order (data points) 
accelerometer_struct.acc_interp_no    = acc_dpint;
accelerometer_struct.sca_interp_no    = sca_dpint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
