function [orbit_model_struct, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori(Xmatrix, sigma_Xgrav, orbit_model_struct, write_gravity_flag) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function: gravity_param_aposteriori.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Gravity Field estimated parameters corrections to background gravity model
%  Write gravity field solution to output data in ICGEM format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - write_gravity_flag  : Gravity sperical harmoncis coefficients to be written
%   1. Estimated Signal (Estimated coefficients corrections)
%   2. Final spherical harmonics coefficients :: Estimated Signal + Background gravity model
% 
% Output arguments: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas D. Papanikolaou                                 25 August 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 14/12/2022, Thomas Loudis Papanikolaou 
%             Code modified to be compatible with structure matrix variables 
% 07/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
% 29/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
% 15/10/2026  Thomas Loudis Papanikolaou 
%             Source Code upgrade  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Old version
% ic_apriori_01 = 1; 
% [orbit_model_struct, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori_0(Xmatrix, ic_apriori_01, orbit_model_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Forces model structure matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
IC_MJD                  = orbit_model_struct.IC_MJD;  
Zo_ICRF_glb             = orbit_model_struct.IC_CRF;  % = [IC_MJDo IC_Zo_vec']; 
Nparam_GLOB             = orbit_model_struct.forces_param_estim_no; 
Nmodel_PARAM_ESTIM_glob = orbit_model_struct.forces_param_estim_yn; 
gravity_struct          = orbit_model_struct.gravity_field; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%-------------------------------------------------------------------------- 
% Xmatrix : Estimated corrections to parameters 
Xmatrix_grav_param = Xmatrix; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Covariance matrix of gravity parameters
%-------------------------------------------------------------------------- 
% % Reference Variance
% % sigma = sqrt(error_matrix' * error_matrix / (n_obs - m_param) );
% sigma0 = 1;
% % Covariance matrix of parameters
% Cx = sigma0^2 * inv_mat(NEQn, inv_id);
% % De-Scaling
% Cx = Kmat * Cx * Kmat;
% 
% sigma_x_square = diag(Cx);
% sigma_x = sqrt(diag(Cx));
% sigma_x_sqrt = sqrt(sigma_x_square);
% sigma_x_grav = sigma_x;

sigma_x_grav = sigma_Xgrav;
%-------------------------------------------------------------------------- 
  
%-------------------------------------------------------------------------- 
% Gravity Field Parameters Corrections to Triangulation matrices
%-------------------------------------------------------------------------- 
C_degree_order = gravity_struct.C_degree_order_estim; 
S_degree_order = gravity_struct.S_degree_order_estim; 
N_param_GRAV = gravity_struct.parameters_number;
[Nparam_C, k] = size(C_degree_order); 
degree_max_estim = C_degree_order(Nparam_C,1) ;

% Gravity Field Parameters 
% [dC_matrix_aposteriori, dS_matrix_aposteriori, Xaposteriori_grav] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, Xmatrix_grav_param); 
% [dC_triang,dS_triang] = gravity_param_vec2triang(C_degree_order, S_degree_order, Xmatrix);
[dC_matrix_aposteriori,dS_matrix_aposteriori] = gravity_param_vec2triang(C_degree_order, S_degree_order, Xmatrix);

% Gravity Field Parameters Covariance matrix
% [sigma_C_matrix, sigma_S_matrix, sigma_Xaposteriori_P] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, sigma_x_grav);
[sigma_C_matrix, sigma_S_matrix] = gravity_param_vec2triang(C_degree_order, S_degree_order, sigma_x_grav);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Gravity field background model
Cnm_ref_model = gravity_struct.Cnm ; 
Snm_ref_model = gravity_struct.Snm ; 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Estimated Signal :: Apriori values (zero or values from previous iteration of the estimator)
delta_Cnm_estim = gravity_struct.Cnm_estim;  
delta_Snm_estim = gravity_struct.Snm_estim;  
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Reference background gravity model to Apriori values
Cnm_apriori = delta_Cnm_estim ; 
Snm_apriori = delta_Snm_estim ; 
%-------------------------------------------------------------------------- 

% grav_background_stat_signal = 1;
% if grav_background_stat_signal == 1
% Update Cnm, Snm matrices :: Add harmonics coefficients corrections to background gravity model
[Cmatrix_aposteriori, Smatrix_aposteriori] = harmonics_sum(Cnm_apriori,Snm_apriori, dC_matrix_aposteriori,dS_matrix_aposteriori,-1);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Update gravity_struct 
% Estimated corrections to gravity parameters
gravity_struct.Cnm_estim  = Cmatrix_aposteriori;  
gravity_struct.Snm_estim  = Smatrix_aposteriori;  

% Gravity parameters (including estimated corrections to coefficients)
% gravity_struct.Cnm = Cmatrix_aposteriori; 
% gravity_struct.Snm = Smatrix_aposteriori; 

% Update Orbit model matrix :: orbit_model_struct
orbit_model_struct.gravity_field = gravity_struct; 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% Overall Gravity signal adding the estimated signal
% (truncated to the maximum degree of the estimated coefficients)
%-------------------------------------------------------------------------- 
% Gravity field background model and estimated gravity signal 
[Cmatrix_sum_signal, Smatrix_sum_signal] = harmonics_sum(Cmatrix_aposteriori, Smatrix_aposteriori, Cnm_ref_model, Cnm_ref_model, degree_max_estim);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Case of (temporal) Gravity signal included in the background; special cases only 
%-------------------------------------------------------------------------- 
% Background Gravity signal (temporal) included in the initial data; Optional, applied only for case studies
gravity_signal_yn = gravity_struct.gravity_signal_yn;

test_gravity_signal_yn    = strcmp(gravity_signal_yn,'y'); 
if test_gravity_signal_yn == 1 

% Signal' delta C,S coefficients matrices
Cnm_signal = gravity_struct.Cnm_signal;
Snm_signal = gravity_struct.Snm_signal;
degree_signal = gravity_struct.degree_signal;
% degree_signal :: [degree_signal order_signal]

% Add temporal gravity signal to harmonics coefficients for the gravity coefficients to output data
[Cmatrix_sum_signal, Smatrix_sum_signal] = harmonics_sum(Cmatrix_sum_signal, Smatrix_sum_signal, Cnm_signal, Snm_signal,degree_max_estim);
end
%-------------------------------------------------------------------------- 



%-------------------------------------------------------------------------- 
% Gravity sperical harmoncis coefficients to be written
% 1. Estimated Signal (Estimated coefficients corrections)
% 2. Final spherical harmonics coefficients :: Estimated Signal + Background gravity model
%-------------------------------------------------------------------------- 
if write_gravity_flag == 1
Cmatrix_model = dC_matrix_aposteriori;
Smatrix_model = dS_matrix_aposteriori;
elseif write_gravity_flag == 2
Cmatrix_model = Cmatrix_sum_signal;
Smatrix_model = Smatrix_sum_signal;
end
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Write Gravity field solution to gravity model icgem format
%-------------------------------------------------------------------------- 
GM = gravity_struct.GM; 
radius = gravity_struct.radius; 

[n,m] = size(Cmatrix_model);
sigma_Cnm = sigma_C_matrix; 
sigma_Snm = sigma_S_matrix;

tide_system = 'tide_free';
time_period_of_data = 'daily';
n_trunc = -1;
gravity_model_name_0 = 'GEORB_gravity_solution';
% gravity_model_filename = 'DORUS_grav_sol'

% gravity_model_date = sprintf('%s%s%s',mjd_start,'-',mjd_end); 
gravity_model_date = sprintf('%d',fix(IC_MJD)); 

gravity_model_ext = '.gfc';

gravity_model_filename = sprintf('%s%s%s%s',gravity_model_name_0,'_',gravity_model_date,gravity_model_ext); 

% [gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cnm,Snm, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)
[gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM, radius, Cmatrix_model, Smatrix_model, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data);
%-------------------------------------------------------------------------