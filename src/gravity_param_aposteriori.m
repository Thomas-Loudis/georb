function [orbit_model_struct, Cmatrix_aposteriori, Smatrix_aposteriori, gravity_model_filename] = gravity_param_aposteriori(Xmatrix, ic_apriori_01, orbit_model_struct)


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
% Dr. Thomas D. Papanikolaou                                 25 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 14/12/2022, Thomas Loudis Papanikolaou
%             Code modified to be compatible with structure matrix variables
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
% 29/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forces model structure matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC_MJD                  = orbit_model_struct.IC_MJD; 
Zo_ICRF_glb             = orbit_model_struct.IC_CRF;  % = [IC_MJDo IC_Zo_vec'];
Nparam_GLOB             = orbit_model_struct.forces_param_estim_no;
Nmodel_PARAM_ESTIM_glob = orbit_model_struct.forces_param_estim_yn;
gravity_struct          = orbit_model_struct.gravity_field;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xmatrix : Estimated corrections to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Estimated corrections: Initial State vector
% Xmatrix_Zo = Xmatrix(1:6,1);
% Estimated corrections: Additional unknown parameters (force related)
% Xmatrix_P = Xmatrix(7:end,1);

Xmatrix_grav_param = Xmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aposteriori matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d1 d2] = size(Xmatrix);
% Nparam_Xmatrix_P = d1 - 6;
% Xaposteriori_Z = zeros(6,1);
% Xaposteriori_P = zeros(Nparam_Xmatrix_P,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Estimated Parameters' set to aposteriori estimated values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown parameters: Force related parameters
if Nparam_GLOB > 0
% Nparam = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Nmodel_PARAM_ESTIM_glob(4) == 1    
if Nmodel_PARAM_ESTIM_glob(4) == 1  && ic_apriori_01 > 0

meth_CS_cor = 1;

% Gravity Field Parameters matrix
if meth_CS_cor == 1
dC_matrix_apriori = gravity_struct.Cnm_estim;
dS_matrix_apriori = gravity_struct.Snm_estim;
elseif meth_CS_cor == 2
dC_matrix_apriori = gravity_struct.Cnm;
dS_matrix_apriori = gravity_struct.Snm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_degree_order = gravity_struct.C_degree_order_estim;
S_degree_order = gravity_struct.S_degree_order_estim;
N_param_GRAV = gravity_struct.parameters_number;
N_param_GRAV_function_gravity_param_aposteriori = gravity_struct.parameters_number

% Xmatrix_grav_param = Xmatrix_P(Nparam + 1 : Nparam + N_param_GRAV);

[dC_matrix_aposteriori, dS_matrix_aposteriori, Xaposteriori_grav] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, Xmatrix_grav_param);
% Xaposteriori_P(Nparam + 1 : Nparam + N_param_GRAV,1) = Xaposteriori_grav(:,1);
% Nparam = Nparam + N_param_GRAV;
delta_X_gravparam = Xaposteriori_grav - Xmatrix_grav_param;
delta_X_gravparam_sum = sum(delta_X_gravparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if meth_CS_cor == 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Cmatrix_aposteriori = dC_matrix_aposteriori;
Smatrix_aposteriori = dS_matrix_aposteriori;

gravity_struct.Cnm = Cmatrix_aposteriori;
gravity_struct.Snm = Smatrix_aposteriori;
% Update Forces model matrix 
orbit_model_struct.gravity_field = gravity_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif meth_CS_cor == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update parameters
approach_no = 1;
% 1. Update Gravity temporal signal delta_Cnm, delta_Snm matrices
if approach_no == 1
gravity_struct.Cnm_estim  = dC_matrix_aposteriori; 
gravity_struct.Snm_estim  = dS_matrix_aposteriori; 
% Update Forces model matrix 
orbit_model_struct.gravity_field = gravity_struct;

Cmatrix_aposteriori = gravity_struct.Cnm_estim;
Smatrix_aposteriori = gravity_struct.Snm_estim;

% 2. temp approach (update of the gravity field model matrices)
elseif approach_no == 2
Cnm_apriori = gravity_struct.Cnm ;
Snm_apriori = gravity_struct.Snm ;
[Cmatrix_aposteriori, Smatrix_aposteriori] = harmonics_sum(Cnm_apriori,Snm_apriori, dC_matrix_aposteriori,dS_matrix_aposteriori,-1);
gravity_struct.Cnm = Cmatrix_aposteriori;
gravity_struct.Snm = Smatrix_aposteriori;
% Update Forces model matrix 
orbit_model_struct.gravity_field = gravity_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


IC_MJD_gravity_param_aposteriori = IC_MJD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Gravity field solution to gravity model icgem format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM = gravity_struct.GM;
radius = gravity_struct.radius;
% TEMP
[n,m] = size(Cmatrix_aposteriori);
sigma_Cnm = zeros(n,n); 
sigma_Snm = sigma_Cnm;

tide_system = 'tide_free';
time_period_of_data = 'daily';
n_trunc = -1;
gravity_model_name_0 = 'GEORB_grav_sol';
gravity_model_date = sprintf('%d',fix(IC_MJD));
gravity_model_ext = '.gfc';
%gravity_model_filename = 'DORUS_GRACE-FO_2021-07-17.gfc'
gravity_model_filename = sprintf('%s%s%s%s',gravity_model_name_0,'_',gravity_model_date,gravity_model_ext);

% [gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cnm,Snm, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)
[gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cmatrix_aposteriori,Smatrix_aposteriori, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ic_apriori_01 == 0
% TEMP     
dC_matrix_aposteriori = gravity_struct.Cnm_estim; 
dS_matrix_aposteriori = gravity_struct.Snm_estim; 

C_matrix_aposteriori = gravity_struct.Cnm_estim; 
S_matrix_aposteriori = gravity_struct.Snm_estim; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
else    
% Aposteriori matrix
% Xaposteriori = Xaposteriori_Z;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT arguments
% Update Forces model matrix 
% orbit_model_struct.gravity_field = gravity_struct;

% Cmatrix_aposteriori = gravity_struct.Cnm_estim;
% Smatrix_aposteriori = gravity_struct.Snm_estim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not supported by current Version 
% ...

% if Nparam_common > 0
if 1 < 0

clear Xmatrix_P Xaposteriori_P   
Xmatrix_P = Xmatrix_gravparam;
% if comb_solution == 4
%     Xmatrix_P = Xmatrix_common;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field Parameters matrix
gravity_struct = orbit_model_matrix_GRACE1.gravity_field;
dC_matrix_apriori = gravity_struct.Cnm_estim;
dS_matrix_apriori = gravity_struct.Snm_estim;
C_degree_order = gravity_struct.C_degree_order_estim;
S_degree_order = gravity_struct.S_degree_order_estim;
% N_param_GRAV = gravity_struct.parameters_number;
% Xmatrix_grav_param = Xmatrix_P(Nparam + 1 : Nparam + N_param_GRAV);
[dC_matrix_aposteriori, dS_matrix_aposteriori] = gravity_param_cor(dC_matrix_apriori, dS_matrix_apriori, C_degree_order, S_degree_order, Xmatrix_P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Cnm, Snm matrices
Cnm_apriori = gravity_struct.Cnm ;
Snm_apriori = gravity_struct.Snm ;
[Cmatrix_aposteriori, Smatrix_aposteriori] = harmonics_sum(Cnm_apriori,Snm_apriori, dC_matrix_aposteriori,dS_matrix_aposteriori,-1);

% gravity_struct.Cnm = Cmatrix_aposteriori;
% gravity_struct.Snm = Smatrix_aposteriori;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Gravity field solution to gravity model icgem format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM = gravity_struct.GM;
radius = gravity_struct.radius;
% TEMP
[n,m] = size(Cmatrix_aposteriori);
sigma_Cnm = zeros(n,n); 
sigma_Snm = sigma_Cnm;

tide_system = 'tide_free';
time_period_of_data = 'daily';
n_trunc = -1;
gravity_model_name_0 = 'DORUS_GRACE-FO';
% MJD_0 = fix(orbcGA(1,1));
% MJD_0 = fix(orbit_model_matrix_GRACE1.IC_MJD);
gravity_model_date = sprintf('%s',IC_MJD);
gravity_model_ext = '.gfc';
%gravity_model_filename = 'DORUS_GRACE-FO_2021-07-17.gfc'
gravity_model_filename = sprintf('%s%s%s%s',gravity_model_name_0,'_',gravity_model_date,gravity_model_ext); 

% [gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cnm,Snm, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)
[gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cmatrix_aposteriori,Smatrix_aposteriori, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_fname = sprintf('%s%d%s','dC_matrix_aposteriori_Nest',i_iter_estim,'.est')
% save out_fname dC_matrix_aposteriori -ASCII -double
save(out_fname, 'dC_matrix_aposteriori', '-ascii')

out_fname = sprintf('%s%d%s','dS_matrix_aposteriori_Nest',i_iter_estim,'.est')
% save out_fname dS_matrix_aposteriori -ASCII -double
save(out_fname, 'dS_matrix_aposteriori', '-ascii')

save dC_matrix.gfc dC_matrix_aposteriori -ASCII -double
save dS_matrix.gfc dS_matrix_aposteriori -ASCII -double

save Cmatrix_aposteriori.gfc Cmatrix_aposteriori -ASCII -double
save Smatrix_aposteriori.gfc Smatrix_aposteriori -ASCII -double

save Cmatrix_aposteriori.out Cmatrix_aposteriori -ASCII -double
save Smatrix_aposteriori.out Smatrix_aposteriori -ASCII -double

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%