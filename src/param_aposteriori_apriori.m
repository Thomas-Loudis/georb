function [Zo_estim, Xaposteriori] = param_aposteriori_apriori(Xmatrix, ic_apriori_01)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: param_aposteriori.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Estimated Parameters: set aposteriori estimate values to parameters' variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - ic_apriori_01 :
%   0 : Apriori values of IC Parameters are provided by the variable Xmatrix
%   1 : Corrections to IC Parameters are provided by the variable Xmatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas D. Papanikolaou                                    August 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 21/04/2021, Dr. Thomas Papanikolaou
%             Function based on code extracted from old mainf_DOD.m function 
% 25/08/2022, Thomas Loudis Papanikolaou
%             Modified to support apriori and aposteriori values and
%             corrections
% 14/12/2022, Thomas Loudis Papanikolaou
%             Code modified to be compatible the global structure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global Zo_ICRF_glb
global Nparam_GLOB 
global Nmodel_PARAM_ESTIM_glob
global emp_cpr_glob
global accelerometer_data_cal_glob
global pulses_stoch_accel_glob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xmatrix : Estimated corrections to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated corrections: Initial State vector
Xmatrix_Zo = Xmatrix(1:6,1);
% Estimated corrections: Additional unknown parameters (force related)
Xmatrix_P = Xmatrix(7:end,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aposteriori matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d1 d2] = size(Xmatrix);
Nparam_Xmatrix_P = d1 - 6;
Xaposteriori_Z = zeros(6,1);
Xaposteriori_P = zeros(Nparam_Xmatrix_P,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Estimated Parameters' set to aposteriori estimated values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Initial State Vector: set a-posteriori estimated values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
Zo_apriori = Zo_ICRF_glb(1,2:7)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic_apriori_01 == 0
    Zo_apriori = zeros(6,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zo_aposteriori = Zo_apriori + Xmatrix_Zo;
Xaposteriori_Z = Zo_aposteriori;
Zo_estim = [Zo_ICRF_glb(1,1) Zo_aposteriori'];
Zo_ICRF_glb = Zo_estim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown parameters: Force related parameters
if Nparam_GLOB > 0
Nparam = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces Cycle-Per-Revolution parameters: Number of unknown parameters to be estimated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias and cycle-per-revolution coefficients for sin & cos terms
if Nmodel_PARAM_ESTIM_glob(1) == 1        
% Obtain variables from global structure array
emp_cpr_struct = emp_cpr_glob;
% Empirical accelerations vector matrix 3x3 for radial, along-track, cross-track
EMP_ACCEL_apriori = emp_cpr_struct.acceleration_matrix; 
[n1 n2] = size(EMP_ACCEL_apriori);
EMP_ACCEL_aposteriori = zeros(n1,n2);
% Array for defining the parameters to be estimated
EMP_FORCE_PARAM_yn = emp_cpr_struct.parameters_01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic_apriori_01 == 0
    EMP_ACCEL_apriori = zeros(n1,n2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if EMP_FORCE_PARAM_yn(1,1) == 1
    Nparam = Nparam + 1;
    % bias_radial = EMP_ACCEL_apriori(1,1) + Xmatrix_P(Nparam);
    EMP_ACCEL_aposteriori(1,1) = EMP_ACCEL_apriori(1,1) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(1,1);
end 
if EMP_FORCE_PARAM_yn(2,1) == 1
    Nparam = Nparam + 1;
    %bias_along 
    EMP_ACCEL_aposteriori(2,1) = EMP_ACCEL_apriori(2,1) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(2,1);
end
if EMP_FORCE_PARAM_yn(3,1) == 1
    Nparam = Nparam + 1;
    % bias_cross 
    EMP_ACCEL_aposteriori(3,1) = EMP_ACCEL_apriori(3,1) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(3,1);
end

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle-per-Revolution (CPR) terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPR radial
if EMP_FORCE_PARAM_yn(1,2) == 1
    Nparam = Nparam + 1;
    % Cr 
    EMP_ACCEL_aposteriori(1,2) = EMP_ACCEL_apriori(1,2) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(1,2);
end

if EMP_FORCE_PARAM_yn(1,3) == 1
    Nparam = Nparam + 1;
    % Sr 
    EMP_ACCEL_aposteriori(1,3) = EMP_ACCEL_apriori(1,3) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(1,3);
end 

% CPR along-track
if EMP_FORCE_PARAM_yn(2,2) == 1
    Nparam = Nparam + 1;
    % Ct 
    EMP_ACCEL_aposteriori(2,2) = EMP_ACCEL_apriori(2,2) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(2,2);
end

if EMP_FORCE_PARAM_yn(2,3) == 1
    Nparam = Nparam + 1;
    % St 
    EMP_ACCEL_aposteriori(2,3) = EMP_ACCEL_apriori(2,3) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(2,3);
end

% CPR cross-track
if EMP_FORCE_PARAM_yn(3,2) == 1
    Nparam = Nparam + 1;
    % Cn 
    EMP_ACCEL_aposteriori(3,2) = EMP_ACCEL_apriori(3,2) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(3,2);
end

if EMP_FORCE_PARAM_yn(3,3) == 1
    Nparam = Nparam + 1;
    % Sn 
    EMP_ACCEL_aposteriori(3,3) = EMP_ACCEL_apriori(3,3) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = EMP_ACCEL_aposteriori(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Parameters values Update :: Empirical Forces global variable 
% EMP_ACCEL_glob = EMP_ACCEL_aposteriori;
emp_cpr_glob.acceleration_matrix = EMP_ACCEL_aposteriori;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE Accelerometry data Calibration Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nmodel_PARAM_ESTIM_glob(2) == 1    
% Accelerometer Calibration Parameters matrix
ACC_CAL_PARAM_matrix = accelerometer_data_cal_glob.cal_parameters;
acc_cal_param_apriori = ACC_CAL_PARAM_matrix;
[d1, d2] = size(ACC_CAL_PARAM_matrix);
acc_cal_param_aposteriori = zeros(d1,d2);
Nparam_ACC_CAL = d1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ic_apriori_01 == 0
    acc_cal_param_apriori = zeros(d1,d2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_cal_param = 1 : Nparam_ACC_CAL
    Nparam = Nparam + 1;
    acc_cal_param_aposteriori(i_cal_param,1) = acc_cal_param_apriori(i_cal_param,1) + Xmatrix_P(Nparam);
    Xaposteriori_P(Nparam,1) = acc_cal_param_aposteriori(i_cal_param,1);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Parameter value update :: CAL parameters global variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACC_CAL_PARAM_matrix = acc_cal_param_aposteriori;
accelerometer_data_cal_glob.cal_parameters = ACC_CAL_PARAM_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Accelerations/Pulses parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nmodel_PARAM_ESTIM_glob(3) == 1    
    pulses_accel_struct = pulses_stoch_accel_glob;
    % Empirical accelerations matrix 
    pulses_matrix = pulses_accel_struct.acceleration_matrix;
    [n1, n2] = size(pulses_matrix);    
    % Apriori parameters' values    
    Nparam_pulses_k = n1 * (n2 - 2);
    pulses_parameters_matrix = zeros(Nparam_pulses_k , 1);
    k = 0;
    for i = 1 : n1
        for j = 3 : n2
            k = k + 1;
            pulses_parameters_matrix(k,1) = pulses_matrix(i,j);
        end
    end
    pulses_param_apriori = pulses_parameters_matrix;
    Nparam_pulses = k;
    
    % Aposteriori parameters' values
    pulses_param_aposteriori = zeros(Nparam_pulses,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ic_apriori_01 == 0
        pulses_param_apriori = pulses_param_aposteriori;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i_param = 1 : Nparam_pulses
        Nparam = Nparam + 1;
        %Xcorr_pulse = Xmatrix_P(Nparam)
        pulses_param_aposteriori(i_param,1) = pulses_param_apriori(i_param,1) + Xmatrix_P(Nparam);
        Xaposteriori_P(Nparam,1) = pulses_param_aposteriori(i_param,1);            
    end    
    
    % Update global Pulses matrix with aposteriori values 
    k = 0;
    for i = 1 : n1
        for j = 3 : n2
            k = k + 1;
            pulses_matrix(i,j) = pulses_param_aposteriori(k,1);  
        end
    end
    pulses_stoch_accel_glob.acceleration_matrix = pulses_matrix;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aposteriori matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xaposteriori = [Xaposteriori_Z ; Xaposteriori_P];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aposteriori matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xaposteriori = [Xaposteriori_Z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
