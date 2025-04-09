function [acc_emp,PD_accemp_Z,PD_accemp_P] = pdv_acclempirical(mjd,rGCRS,vGCRS, nCPR, GM, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
% Empirical accelerations modelling based on CPR terms
% Partial derivatives of empirical forces w.r.t. position/velocity vector 
% and the parameters to be estimated
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - tmjd  :
% - rGCRS : position vector in the inertial frame r_crf = [x y z]'
% - vGCRS : velocity vector in the inertial frame v_crf = [Vx Vy Vz]'
% 
% Output arguments:
% - PD_acc_P : Partial derivatives matrix 3xNp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas D. Papanikolaou                                    5 July 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 14/12/2022, Dr. Thomas Loudis Papanikolaou
%             Code upgrade to use the orbit model structure arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Forces model matrix 
emp_cpr_glob = orbit_model_struct.empirical_forces_cpr; 
accelerometer_data_cal_glob = orbit_model_struct.accelerometer_struct; 
MJDo_glb = orbit_model_struct.IC_MJD;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables of Empirical Forces based on CPR terms
emp_cpr_struct = emp_cpr_glob;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Frame options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Orbital Frame : orbital_frame 
% - Spacecraft reference frame : srf_frame
% reference_frame = EMP_REFRAME_glob;
reference_frame = emp_cpr_struct.reference_frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical forces :: CPR modelling to be considered (y/n effects)
% Array for defining the parameters to be estimated
EMP_FORCE_PARAM_yn = emp_cpr_struct.parameters_01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias accelerations and CPR coefficients for the 3 directions (radial,
% along-track, cross-track) 
% EMP_Accel_matrix = EMP_ACCEL_glob;
EMP_Accel_matrix = emp_cpr_struct.acceleration_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nparam_EMP_FORCE = emp_cpr_struct.parameters_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read bias and cycle-per-revolution coefficients per orbital direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias_radial = 0;
Cr = 0;
Sr = 0;
bias_along = 0;
Ct = 0;
St = 0;
bias_cross = 0;
Cn = 0;
Sn = 0;

if 1<0
    
%if EMP_Force(1,2) == 1
if EMP_FORCE_PARAM_yn(1,1) == 1    
    bias_radial = EMP_Accel_matrix(1,1);
end
if EMP_FORCE_PARAM_yn(1,2) == 1    
    Cr = EMP_Accel_matrix(1,2);
end
if EMP_FORCE_PARAM_yn(1,3) == 1    
    Sr = EMP_Accel_matrix(1,3);
end

%if EMP_Force(1,3) == 1
if EMP_FORCE_PARAM_yn(2,1) == 1    
    bias_along = EMP_Accel_matrix(2,1);
end
if EMP_FORCE_PARAM_yn(2,2) == 1    
    Ct = EMP_Accel_matrix(2,2);
end
if EMP_FORCE_PARAM_yn(3,3) == 1    
    St = EMP_Accel_matrix(2,3);
end

%if EMP_Force(1,4) == 1
if EMP_FORCE_PARAM_yn(3,1) == 1    
    bias_cross = EMP_Accel_matrix(3,1);
end
if EMP_FORCE_PARAM_yn(3,2) == 1    
    Cn = EMP_Accel_matrix(3,2);
end
if EMP_FORCE_PARAM_yn(3,3) == 1    
    Sn = EMP_Accel_matrix(3,3);
end

end

bias_radial = EMP_Accel_matrix(1,1);
Cr          = EMP_Accel_matrix(1,2);
Sr          = EMP_Accel_matrix(1,3);

bias_along  = EMP_Accel_matrix(2,1);
Ct          = EMP_Accel_matrix(2,2);
St          = EMP_Accel_matrix(2,3);
    
bias_cross  = EMP_Accel_matrix(3,1);
Cn          = EMP_Accel_matrix(3,2);
Sn          = EMP_Accel_matrix(3,3);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Frame:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Frame
test = strcmp(reference_frame,'orbital_frame');
if test == 1
% Transformation matrix from inertial to orbital frame
[Rrtn,er,et,en] = orbital_transf(rGCRS,vGCRS);
% Transformation matrix from orbital to inertial frame
Rrtn2crf = inv(Rrtn);
Rmatrix = Rrtn2crf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Spacecraft reference frame (SRF) : GRACE & GRACE-FO missions 
test = strcmp(reference_frame,'srf_frame');
if test == 1    
% sca_array = SCA_array_glob;
% scadpint = SCA_dpint_glob;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation matrix from SRF (Spacecraft reference frame) to ICRF (Celestial
% Reference Frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation matrix based on Quaternion
%q = sca_interp_matrix(i,3:6)';
quat = quaternions_int;
[R_quat] = quat_rot(quat);
% SRF to Celestial reference frame tranfromation matrix
Rsrf2crf = inv(R_quat);
Rmatrix = Rsrf2crf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_bias_rtn = [bias_radial bias_along bias_cross]';
% Transformation to Inertial system
[a_emp_bias] = Rmatrix * a_bias_rtn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle-per-revolution (CPR) empirical accelerations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time dependent argument of CPR cos and sin terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keplerian elements
[a_semiaxis,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_deg] = kepler(rGCRS,vGCRS,GM);
u_rad = u_deg * (pi / 180);
M_rad = M_deg * (pi / 180);

% Orbital Period
T_orb_period = 2*pi * sqrt(a_semiaxis^3 / GM);

% Frequency
freq_orb = 2*pi / T_orb_period;  

% Time difference since start of the orbit arc
mjd_to = MJDo_glb;
delta_mjd = mjd - mjd_to;
delta_tsec = delta_mjd * 24 * 3600;

% Time argument based on time since start and orbital period
freq_t_product = freq_orb * delta_tsec;

% Time-variable argument of CPR terms
cpr_time_argument = u_rad;
cpr_time_argument = M_rad;
cpr_time_argument = freq_t_product;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-cycle-per-revolution (1-CPR) acclerations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-CPR radial    
a_R_cpr = Cr * cos(nCPR * cpr_time_argument) + Sr * sin(nCPR * cpr_time_argument);

% 1-CPR along-track    
a_T_cpr = Ct * cos(nCPR * cpr_time_argument) + St * sin(nCPR * cpr_time_argument);

% 1-CPR cross-track    
a_N_cpr = Cn * cos(nCPR * cpr_time_argument) + Sn * sin(nCPR * cpr_time_argument);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle-per-revolution acceleration vector: transformation to inertial frame
a_cpr_rtn = [a_R_cpr a_T_cpr a_N_cpr]';
a_emp_cpr = Rmatrix * a_cpr_rtn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives of CPR w.r.t. XYZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD_acpr_xyz = [PD_aCPRx_x PD_aCPRx_y PD_aCPRx_z
%                PD_aCPRy_x PD_aCPRy_y PD_aCPRy_z
%                PD_aCPRz_x PD_aCPRz_y PD_aCPRz_z ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d(u)/d(xyz) derivatives through Differential approximations
dxyz = 0.001;
rGCRS_dx = [rGCRS(1,1)+dxyz; rGCRS(2,1); rGCRS(3,1) ];
rGCRS_dy = [rGCRS(1,1); rGCRS(2,1)+dxyz; rGCRS(3,1) ];
rGCRS_dz = [rGCRS(1,1); rGCRS(2,1); rGCRS(3,1)+dxyz ];

[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dx_deg] = kepler(rGCRS_dx,vGCRS,GM);
[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dy_deg] = kepler(rGCRS_dy,vGCRS,GM);
[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dz_deg] = kepler(rGCRS_dz,vGCRS,GM);
u_dx_rad = u_dx_deg * (pi / 180);
u_dy_rad = u_dy_deg * (pi / 180);
u_dz_rad = u_dz_deg * (pi / 180);

PD_u_x = (u_dx_rad - u_rad) / dxyz;
PD_u_y = (u_dy_rad - u_rad) / dxyz;
PD_u_z = (u_dz_rad - u_rad) / dxyz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Partial derivatives w.r.t. Velocity vector 
dVxyz = 10^-5;
vGCRS_dVx = [vGCRS(1,1)+dVxyz; vGCRS(2,1); vGCRS(3,1) ];
vGCRS_dVy = [vGCRS(1,1); vGCRS(2,1)+dVxyz; vGCRS(3,1) ];
vGCRS_dVz = [vGCRS(1,1); vGCRS(2,1); vGCRS(3,1)+dVxyz ];

[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dVx_deg] = kepler(rGCRS,vGCRS_dVx,GM);
[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dVy_deg] = kepler(rGCRS,vGCRS_dVy,GM);
[a_deg,e_deg,i_deg,Omega_deg,omega_deg,f_deg,M_deg,E_deg,u_dVz_deg] = kepler(rGCRS,vGCRS_dVz,GM);
u_dVx_rad = u_dVx_deg * (pi / 180);
u_dVy_rad = u_dVy_deg * (pi / 180);
u_dVz_rad = u_dVz_deg * (pi / 180);

PD_u_Vx = (u_dVx_rad - u_rad) / dVxyz;
PD_u_Vy = (u_dVy_rad - u_rad) / dVxyz;
PD_u_Vz = (u_dVz_rad - u_rad) / dVxyz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%PD_aCPRx_u = (-et(1,1) * Ct * sin(nCPR * u_rad) * nCPR + et(1,1) * St * cos(nCPR * u_rad)*nCPR);
%PD_aCPRy_u = (-et(2,1) * Ct * sin(nCPR * u_rad) * nCPR + et(2,1) * St * cos(nCPR * u_rad)*nCPR);
%PD_aCPRz_u = (-et(3,1) * Ct * sin(nCPR*u_rad)*nCPR + et(3,1) * St * cos(nCPR*u_rad)*nCPR);

%PD_acpr_xyz = [ (PD_aCPRx_u * PD_u_x)  (PD_aCPRx_u * PD_u_y)  (PD_aCPRx_u * PD_u_z)
%                (PD_aCPRy_u * PD_u_x)  (PD_aCPRy_u * PD_u_y)  (PD_aCPRy_u * PD_u_z)
%                (PD_aCPRz_u * PD_u_x)  (PD_aCPRz_u * PD_u_y)  (PD_aCPRz_u * PD_u_z) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives: d(acpr)icrf/d(xyz) = Rrtn' * d(acpr)rtn/d(xyz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%d(acpr)rtn/d(u)
pd_aR_u = -Cr * sin(nCPR * cpr_time_argument) * nCPR + Sr * cos(nCPR * cpr_time_argument) * nCPR; 
pd_aT_u = -Ct * sin(nCPR * cpr_time_argument) * nCPR + St * cos(nCPR * cpr_time_argument) * nCPR;
pd_aN_u = -Cn * sin(nCPR * cpr_time_argument) * nCPR + Sn * cos(nCPR * cpr_time_argument) * nCPR;

%d(acpr)rtn/d(xyz)
PD_acpr_RTN_xyz(1,1:3) = pd_aR_u * [PD_u_x PD_u_y PD_u_z];
PD_acpr_RTN_xyz(2,1:3) = pd_aT_u * [PD_u_x PD_u_y PD_u_z];
PD_acpr_RTN_xyz(3,1:3) = pd_aN_u * [PD_u_x PD_u_y PD_u_z];
% Transformation to inertial frame
PD_acpr_xyz = Rmatrix * PD_acpr_RTN_xyz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical accelerations sum
acc_emp = a_emp_bias + a_emp_cpr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical accelerations partial derivatives w.r.t. position vector 
%PD_accemp_Z = PD_acpr_xyz;
PD_accemp_Z = zeros(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical accelerations partial derivatives w.r.t. unknow parameters to be estimated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD_accemp_P = zeros(3,6);    
PD_accemp_P = zeros(3,Nparam_EMP_FORCE);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nparam = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if EMP_FORCE_Bias_yn_glob == 1

%if EMP_Force(1,2) == 1
if EMP_FORCE_PARAM_yn(1,1) == 1    
Nparam = Nparam + 1;
% Bias Radial
pd_aRTN_biasR = [1; 0; 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_biasR;
end

%if EMP_Force(1,3) == 1
if EMP_FORCE_PARAM_yn(2,1) == 1    
Nparam = Nparam + 1;
% Bias Along-track
pd_aRTN_biasT = [0; 1; 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_biasT;
end

%if EMP_Force(1,4) == 1
if EMP_FORCE_PARAM_yn(3,1) == 1    
Nparam = Nparam + 1;
% Bias Cross-track
pd_aRTN_biasN = [0; 0; 1];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_biasN;
end

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPR parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPR Radial
%if EMP_Force(1,2) == 1
if EMP_FORCE_PARAM_yn(1,2) == 1    
Nparam = Nparam + 1;
% Cr radial
pd_aRTN_Cr = [cos(nCPR * cpr_time_argument); 0; 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_Cr;
end

if EMP_FORCE_PARAM_yn(1,3) == 1    
Nparam = Nparam + 1;
% Sr radial
pd_aRTN_Sr = [sin(nCPR * cpr_time_argument); 0; 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_Sr;
end


% CPR Along-track
%if EMP_Force(1,3) == 1
if EMP_FORCE_PARAM_yn(2,2) == 1    
Nparam = Nparam + 1;
% Ct along-track
pd_aRTN_Ct = [0; cos(nCPR * cpr_time_argument); 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_Ct;
end

if EMP_FORCE_PARAM_yn(2,3) == 1    
Nparam = Nparam + 1;
% St along-track
pd_aRTN_St = [0; sin(nCPR * cpr_time_argument); 0];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_St;
end


% CPR Cross-track
%if EMP_Force(1,4) == 1
if EMP_FORCE_PARAM_yn(3,2) == 1    
Nparam = Nparam + 1;
% Cn cross-track
pd_aRTN_Cn = [0; 0; cos(nCPR * cpr_time_argument)];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_Cn;
end

if EMP_FORCE_PARAM_yn(3,3) == 1    
Nparam = Nparam + 1;
% Sn cross-track
pd_aRTN_Sn = [0; 0; sin(nCPR * cpr_time_argument)];
PD_accemp_P(:,Nparam) = Rmatrix * pd_aRTN_Sn;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

