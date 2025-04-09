function [PD_ACC_Cal_Param, facc_ICRF, facc_SRF] = pdv_grace_acc(mjd, acc_cal_param, acc1b_array_TT, sca_array, acc_dpint, scadpint, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: pdv_grace_acc.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Partial derivatives of accelerometry data w.r.t. the calibration parameters
%  to be estimated and the position/velocity vector
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - tmjd  :
% - rGCRS : position vector in the inertial frame r_crf = [x y z]'
% - vGCRS : velocity vector in the inertial frame v_crf = [Vx Vy Vz]'
% 
% Output arguments:
% - PD_ACC_Cal : Partial derivatives matrix 3xNp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Papanikolaou                                     22 April 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 08/07/2022  Thomas Loudis Papanikolaou
%             Code minor modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Forces model matrix 
accelerometer_data_cal_glob = orbit_model_struct.accelerometer_struct; 
MJDo_glb = orbit_model_struct.IC_MJD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accelerometer_struct = accelerometer_data_cal_glob;
acc_cal_bias_yn      = accelerometer_struct.cal_bias_yn;
acc_cal_scale_type   = accelerometer_struct.cal_scale_type;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Calibration Parameters
[d1, d2] = size(acc_cal_param);
Nparam_ACC_CAL_model = d1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration modelling parameters set included:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias parameters
bias_yn         = acc_cal_bias_yn(1);

% Bias Drift (1st order)
bias_drift_1_yn = acc_cal_bias_yn(2);

% Bias Drift (2nd order)
bias_drift_2_yn = acc_cal_bias_yn(3);

% Scale factors matrix
scale_matrix_type = acc_cal_scale_type;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time difference since start of the orbit arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mjd_to = MJDo_glb;
delta_mjd = mjd - mjd_to;
delta_tsec = delta_mjd * 24 * 3600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acc1b Accelerometry data interpolation to the computation epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[acc1b_X_int] = interp_Lagrange(acc1b_array_TT(:,1),acc1b_array_TT(:,3),mjd,acc_dpint);
[acc1b_Y_int] = interp_Lagrange(acc1b_array_TT(:,1),acc1b_array_TT(:,4),mjd,acc_dpint);
[acc1b_Z_int] = interp_Lagrange(acc1b_array_TT(:,1),acc1b_array_TT(:,5),mjd,acc_dpint);
acc1b_int = [acc1b_X_int; acc1b_Y_int; acc1b_Z_int];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix initialisation
facc_SRF = zeros(3,1);        
Nparam = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Calibration modelling to accelerometry observations (acc1b data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias Terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bias_yn == 1
    facc_SRF = [(acc_cal_param(1,1) ); ... 
                (acc_cal_param(2,1) ); ... 
                (acc_cal_param(3,1) )];
        
    Nparam = Nparam + 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = [ Sx  Sxy Sxz 
%       Syx Sy  Syz
%       Szx Szy Sz  ]

ax1b = acc1b_int(1,1);
ay1b = acc1b_int(2,1);
az1b = acc1b_int(3,1);

SCALE_test = 0;

% Diagonal Scale matrix
SCALE_test = strcmp(scale_matrix_type,'diagonal');
if SCALE_test == 1 
    Nparam = Nparam + 1;
    Sx  = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sy  = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sz  = acc_cal_param(Nparam,1);

    Sxy = 0;
    Sxz = 0;
    Syx = 0;
    Syz = 0;
    Szx = 0;
    Szy = 0;    
end

SCALE_test = strcmp(scale_matrix_type,'semi-full');
if SCALE_test == 1 
    Nparam = Nparam + 1;
    Sx  = acc_cal_param(Nparam,1);

    Nparam = Nparam + 1;
    Sxy = acc_cal_param(Nparam,1);

    Nparam = Nparam + 1;
    Sxz = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sy  = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Syz = acc_cal_param(Nparam,1);

    Nparam = Nparam + 1;
    Sz  = acc_cal_param(Nparam,1);

    Syx = Sxy;
    Szx = Sxz;
    Szy = Syz;    
end

% Full Scale matrix 3x3 
SCALE_test = strcmp(scale_matrix_type,'full');
if SCALE_test == 1
%if scale_full_matrix_yn == 1
    Nparam = Nparam + 1;
    Sx  = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sxy = acc_cal_param(Nparam,1);

    Nparam = Nparam + 1;
    Sxz = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Syx = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sy  = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Syz = acc_cal_param(Nparam,1);

    Nparam = Nparam + 1;
    Szx = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Szy = acc_cal_param(Nparam,1);
    
    Nparam = Nparam + 1;
    Sz  = acc_cal_param(Nparam,1);    
end

% Scale terms acceleration in SRF 
facc_scale_full_matrix(1,1) = [Sx * ax1b  + Sxy * ay1b + Sxz * az1b];
facc_scale_full_matrix(2,1) = [Syx * ax1b + Sy * ay1b  + Syz * az1b];
facc_scale_full_matrix(3,1) = [Szx * ax1b + Szy * ay1b + Sz * az1b ];
facc_SRF = facc_SRF + facc_scale_full_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias drift 1st order        
if bias_drift_1_yn == 1
    bias_drift_1_matrix(1,1) = [acc_cal_param(Nparam+1,1) * delta_tsec];
    bias_drift_1_matrix(2,1) = [acc_cal_param(Nparam+2,1) * delta_tsec];
    bias_drift_1_matrix(3,1) = [acc_cal_param(Nparam+3,1) * delta_tsec];
    facc_SRF = facc_SRF + bias_drift_1_matrix;        
    Nparam = Nparam + 3;
end

% bias drift 2nd order        
if bias_drift_2_yn == 1
    bias_drift_2_matrix(1,1) = [acc_cal_param(Nparam+1,1) * delta_tsec^2];
    bias_drift_2_matrix(2,1) = [acc_cal_param(Nparam+2,1) * delta_tsec^2];
    bias_drift_2_matrix(3,1) = [acc_cal_param(Nparam+3,1) * delta_tsec^2];
    facc_SRF = facc_SRF + bias_drift_2_matrix;        
    Nparam = Nparam + 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
% Accleration vector Transformation SRF to ICRF
facc_icrs = Rsrf2crf * facc_SRF;
facc_ICRF = facc_icrs;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometry data partial derivatives w.r.t. calibration parameters to
% be estimated e.g. biases and scale factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
Nparam_Acc_Cal = Nparam_ACC_CAL_model;
PD_ACC_Cal_Param = zeros(3,Nparam_Acc_Cal);    

% Partials in the SRF (spacecraft reference frame)
PD_ACC_Cal_SRF  = zeros(3,Nparam_Acc_Cal);

% Partials in ICRF (Celestial reference frame)
PD_ACC_Cal_ICRF = zeros(3,Nparam_Acc_Cal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nparam = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partials w.r.t. Bias parameters per direction XYZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias acceleration terms        
if bias_yn == 1
    % bias_x
    PD_ACC_Cal_SRF(:,Nparam+1) = [1; 0; 0];
    % bias_y
    PD_ACC_Cal_SRF(:,Nparam+2) = [0; 1; 0];
    % bias_z
    PD_ACC_Cal_SRF(:,Nparam+3) = [0; 0; 1];
    
    Nparam = Nparam + 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partials w.r.t. Scale factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCALE_test = 0;

% Diagonal Scale matrix
SCALE_test = strcmp(scale_matrix_type,'diagonal');
if SCALE_test == 1 
    % Sx
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [ax1b; 0; 0];

    % Sy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; ay1b; 0];

    % Sz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; 0; az1b];       
end


% Semi-Full Scale matrix [Sxy=Syx Sxz=Szx Syz=Szy]
SCALE_test = strcmp(scale_matrix_type,'semi-full');
if SCALE_test == 1 
    % Sx
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [ax1b; 0; 0];

    % Sxy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [ay1b; ax1b; 0];
    
    % Sxz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [az1b; 0; ax1b];

    % Sy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; ay1b; 0];

    % Syz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; az1b; ay1b];

    % Sz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; 0; az1b];       
end


% Full Scale matrix 3x3 
SCALE_test = strcmp(scale_matrix_type,'full');
if SCALE_test == 1
%if scale_full_matrix_yn == 1                
    % Sx
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [ax1b; 0; 0];

    % Sxy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [ay1b; 0; 0];

    % Sxz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [az1b; 0; 0];

    
    % Syx
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; ax1b; 0];

    % Sy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; ay1b; 0];

    % Syz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; az1b; 0];

    
    % Szx
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; 0; ax1b];
    
    % Szy
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; 0; ay1b];

    % Sz
    Nparam = Nparam + 1;
    PD_ACC_Cal_SRF(:,Nparam) = [0; 0; az1b];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partials w.r.t. Drift of the bias parameters per direction XYZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias drift 1st order        
if bias_drift_1_yn == 1
    % drift_x
    PD_ACC_Cal_SRF(:,Nparam+1) = [delta_tsec; 0; 0];
    % drift_y
    PD_ACC_Cal_SRF(:,Nparam+2) = [0; delta_tsec; 0];
    % drift_z
    PD_ACC_Cal_SRF(:,Nparam+3) = [0; 0; delta_tsec];
    
    Nparam = Nparam + 3;
end

% bias drift 2nd order        
if bias_drift_2_yn == 1
    % drift_x
    PD_ACC_Cal_SRF(:,Nparam+1) = [delta_tsec^2; 0; 0];
    % drift_y
    PD_ACC_Cal_SRF(:,Nparam+2) = [0; delta_tsec^2; 0];
    % drift_z
    PD_ACC_Cal_SRF(:,Nparam+3) = [0; 0; delta_tsec^2];
    
    Nparam = Nparam + 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partials matrix rotation from SRF to ICRF
PD_ACC_Cal_ICRF = Rsrf2crf * PD_ACC_Cal_SRF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final partials matrix
PD_ACC_Cal_Param = PD_ACC_Cal_ICRF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
