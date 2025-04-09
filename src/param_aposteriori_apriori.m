function [Zo_estim, Xaposteriori, orbit_model_struct] = param_aposteriori_apriori(Xmatrix, ic_apriori_01, orbit_model_struct)


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
%             Code modified to be compatible with structure matrix variables
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forces model structure matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC_MJDo = orbit_model_struct.IC_MJD; % = IC_MJDo;
Zo_ICRF_glb = orbit_model_struct.IC_CRF; % = [IC_MJDo IC_Zo_vec'];
Nparam_GLOB = orbit_model_struct.forces_param_estim_no;
Nmodel_PARAM_ESTIM_glob = orbit_model_struct.forces_param_estim_yn;
emp_cpr_glob = orbit_model_struct.empirical_forces_cpr; 
accelerometer_data_cal_glob = orbit_model_struct.accelerometer_struct; 
pulses_stoch_accel_glob = orbit_model_struct.empirical_forces_pulses;
gfm_struct_glob = orbit_model_struct.gravity_field;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Update Forces model matrix 
orbit_model_struct.IC_CRF = Zo_ICRF_glb;
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
% Obtain variables via structure array
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
% Estimated Parameters values Update :: Empirical Forces structure matrix  
% EMP_ACCEL_glob = EMP_ACCEL_aposteriori;
emp_cpr_glob.acceleration_matrix = EMP_ACCEL_aposteriori;
% Update Forces model matrix 
orbit_model_struct.empirical_forces_cpr = emp_cpr_glob; 
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
% Estimated Parameter value update :: CAL parameters via central structure matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACC_CAL_PARAM_matrix = acc_cal_param_aposteriori;
accelerometer_data_cal_glob.cal_parameters = ACC_CAL_PARAM_matrix;
% Update Forces model matrix 
orbit_model_struct.accelerometer_struct = accelerometer_data_cal_glob; 
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
    
    % Update Pulses matrix with aposteriori values 
    k = 0;
    for i = 1 : n1
        for j = 3 : n2
            k = k + 1;
            pulses_matrix(i,j) = pulses_param_aposteriori(k,1);  
        end
    end
    pulses_stoch_accel_glob.acceleration_matrix = pulses_matrix;
    % Update Forces model matrix 
    orbit_model_struct.empirical_forces_pulses = pulses_stoch_accel_glob;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nmodel_PARAM_ESTIM_glob(4) == 1    
% if Nmodel_PARAM_ESTIM_glob(4) == 1  && ic_apriori_01 > 0
if ic_apriori_01 > 0
    
% Gravity Field Parameters matrix
dC_matrix_apriori = gfm_struct_glob.Cnm_estim;
dS_matrix_apriori = gfm_struct_glob.Snm_estim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_degree_order = gfm_struct_glob.C_degree_order_estim;
S_degree_order = gfm_struct_glob.S_degree_order_estim;
N_param_GRAV = gfm_struct_glob.parameters_number;
Xmatrix_grav_param = Xmatrix_P(Nparam + 1 : Nparam + N_param_GRAV);
[dC_matrix_aposteriori, dS_matrix_aposteriori, Xaposteriori_grav] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, Xmatrix_grav_param);
Xaposteriori_P(Nparam + 1 : Nparam + N_param_GRAV,1) = Xaposteriori_grav(:,1);
Nparam = Nparam + N_param_GRAV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update parameters
approach_no = 2;
% 1.
if approach_no == 1
gfm_struct_glob.Cnm_estim  = dC_matrix_aposteriori; 
gfm_struct_glob.Snm_estim  = dS_matrix_aposteriori; 
% 2. temp approach (update of the gravity field model matrices)
elseif approach_no == 2
Cnm_apriori = gfm_struct_glob.Cnm ;
Snm_apriori = gfm_struct_glob.Snm ;
[Cmatrix_aposteriori, Smatrix_aposteriori] = harmonics_sum(Cnm_apriori,Snm_apriori, dC_matrix_aposteriori,dS_matrix_aposteriori,-1);
gfm_struct_glob.Cnm = Cmatrix_aposteriori;
gfm_struct_glob.Snm = Smatrix_aposteriori;
% Update Forces model matrix 
orbit_model_struct.gravity_field = gfm_struct_glob;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif ic_apriori_01 == 0
dC_matrix_aposteriori = gfm_struct_glob.Cnm_estim; 
dS_matrix_aposteriori = gfm_struct_glob.Snm_estim; 

end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aposteriori matrix
Xaposteriori = [Xaposteriori_Z ; Xaposteriori_P];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else    
% Aposteriori matrix
Xaposteriori = Xaposteriori_Z;

end
