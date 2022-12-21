function [acc_cal_param] = grace_acc_cal_init(acc,acc_cal_model,scale_matrix_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_acc_cal_init:  GRACE accelerometry calibration parameters initialisation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Initialisation of the GRACE accelerometry data calibration parameters
%  prior estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - acc:      accelerometry data acc1b
%
% Output arguments:
% - acc_cal_param:   Accelerometry calibration parametes 
%   acc_cal_param = [biasx biasy biasz sx sy sz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 22/04/2021, Dr. Thomas Papanikolaou
%             Code extracted from grace_accproc.m function into an
%             individual function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial Accelerometry data array
% acc rate 1 sec
[sz1 sz2] = size(acc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometry calibration parameters matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerometer Calibration model: effects included parameters 

% Number of calibration parameters
% Bias
Nparam_CAL = 3 * acc_cal_model(1);

% Bias drift matrix
Nparam_CAL = Nparam_CAL + 3 * acc_cal_model(2);

% Bias drift 2nd order matrix
Nparam_CAL = Nparam_CAL + 3 * acc_cal_model(3);

% Scale factors
%scale_matrix_type = SCALE_Matrix_Type_GLOB
SCALE_test = 0;
SCALE_case = 0;
% Diagonal Scale matrix
SCALE_test = strcmp(scale_matrix_type,'diagonal');
if SCALE_test == 1 
    Nparam_CAL = Nparam_CAL + 3;
    SCALE_case = 1;
end
% Semi-Full Scale matrix [Sxy=Syx Sxz=Szx Syz=Szy]
SCALE_test = strcmp(scale_matrix_type,'semi-full');
if SCALE_test == 1 
    Nparam_CAL = Nparam_CAL + 6;
    SCALE_case = 2;
end
% Full Scale matrix 3x3 
SCALE_test = strcmp(scale_matrix_type,'full');
if SCALE_test == 1
    Nparam_CAL = Nparam_CAL + 9;    
    SCALE_case = 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation to zero matrix
acc_cal_param = zeros(Nparam_CAL,1);

Nparam = 0;

% Bias terms
bias_yn = acc_cal_model(1);
if bias_yn == 1
    Nparam = 3 * acc_cal_model(1);
end

% Apriori values
scale_ij_apriori = 10^-5;

% Scale factors initialise to 1
if SCALE_case == 1
    acc_cal_param(Nparam+1,1) = 1;
    acc_cal_param(Nparam+2,1) = 1;
    acc_cal_param(Nparam+3,1) = 1;
    
elseif SCALE_case == 2
    acc_cal_param(Nparam+1 : Nparam+6 , 1) = scale_ij_apriori;
    acc_cal_param(Nparam+1 ,1) = 1;
    acc_cal_param(Nparam+4 ,1) = 1;
    acc_cal_param(Nparam+6 ,1) = 1;    
    
elseif SCALE_case == 3
    acc_cal_param(Nparam+1 : Nparam+9 ,1) = scale_ij_apriori;
    acc_cal_param(Nparam+1 , 1)  = 1;
    acc_cal_param(Nparam+5 ,1)  = 1;
    acc_cal_param(Nparam+9 ,1) = 1;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input (apriori) calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cal_opt = 0;

if cal_opt == 1
    %
    %
elseif cal_opt == 2    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apriori calibration parameters Input apriori calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated parameters from previous dynamic orbit determination
% GRACE-FO 1 (GRACE-C)
if grace == 1
    % Calibration model :: 
    % acc_cal_param_apriori =  
% GRACE-FO 2 (GRACE-D)
elseif grace == 2
    % Calibration model :: 
    % acc_cal_param_apriori =     
end
 
elseif cal_opt == 3    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation to apriori values:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

