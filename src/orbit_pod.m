function [rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,Xmatrix,orbc,orbk,orbt,veqZarray,veqParray,rms_orbc_obs, OBS_matrix, Xaposteriori, OBS_residuals, extorb_crf, forces_accel, Gmatrix, Rmatrix, orbit_model_struct] = orbit_pod (orbit_config_fname, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_pod.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Main function for dynamic orbit modelling and determination.
%  The following oblectives are computed hereby: 
%  - orbit propagation
%  - orbit integration
%  - VEQ (variational equations) integration
%  - VEQ arrays (i.e. State Transition matrix and Sensitivity matrix)
%    computation and storing in .out files
%  - Parameters and Orbit Estimation (DOD algorithm)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in in format 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         July 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 15/04/2021, Thomas Papanikolaou
%             Function renamed & upgraded from mainf20120618.m to DOD.m
% 06/07/2022, Thomas Loudis Papanikolaou
%             Function upgraded and renamed from DOD to orbit_pod
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GM_glob = orbit_model_struct.GM_Earth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit modelling :: Initial Conditions 
% orbit_model (orbit_config_fname);
% [orbit_model_struct] = orbit_model_ic (orbit_config_fname, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination mode
param_keyword = 'orbit_mode';
[orbit_mode] = read_param_cfg(orbit_config_fname,param_keyword);

param_estim_01 = 0;

test_pod_mode = strcmp(orbit_mode,'orbit_determination');
if test_pod_mode == 0    
test = strcmp(orbit_mode,'orbit_propagation_veq');
if test == 1
    %fprintf('\n%s\n','Mode: Orbit Propagation EQM and VEQ');    
    Niter_estim = 0;
end

test = strcmp(orbit_mode,'orbit_propagation_eqm');
if test == 1
    %fprintf('\n%s\n','Mode: Orbit Propagation EQM');    
    Niter_estim = -1;
end
%    fprintf('\n\n%s\n','Mode: Orbit Propagation');
%    Niter_estim = -1;

elseif test_pod_mode == 1
    %fprintf('\n%s\n','Mode: Orbit Determinaton');    
    param_keyword = 'estimator_iterations';
    [estimator_iterations] = read_param_cfg(orbit_config_fname, param_keyword);
    estim_no_iterations = sscanf(estimator_iterations,'%d'); 
    Niter_estim = estim_no_iterations;
    param_estim_01 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-Observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if test_pod_mode == 1
% Pseudo-Observations based on kinematic orbit data
% [obsorbc,obsorbt,obsorbc_ext,obsorbt_ext,obsorbc_full,obsorbt_full,COVobs,COVPform] = orbit_obs(orbit_config_fname);

% observation_struct = observation_data_glob;
% orbit model strucutre matrix
observation_struct = orbit_model_struct.observation_struct;
% Orbit observations in ICRF 
obsorbc = observation_struct.obs_orbit_crf ;
% Orbit observations in ITRF 
obsorbt = observation_struct.obs_orbit_trf ;
% Observations' covariance matrix 
COVobs = observation_struct.obs_cov ;
% Observations' weight matrix form 
COVPform = observation_struct.obs_weight;

% Observations screening based on outlier threshold 
% Read configuration file for pseudo-observations outlier threshold    
param_keyword = 'obs_outliers_yn';
[obs_outliers_yn] = read_param_cfg(orbit_config_fname,param_keyword);

param_keyword = 'pseudo_obs_sigma';
[param_value] = read_param_cfg(orbit_config_fname,param_keyword);
% Pseudo-observations outlier threshold  :: Convert to meters
pseudo_obs_sigma_outlier = sscanf(param_value,'%f %*') / 100;
else
   obsorbc = zeros(4,4); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions: Initial Epoch and State Vector (apriori values)
[orbit_model_struct] = ic_apriori (orbit_config_fname, obsorbc, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Orbit Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination: Orbital Equations numerical integration & Parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test_pod_mode == 1
    %fprintf('\n%s \n\n', 'Dynamic Orbit Determination: in-progress : Iterations of Orbit Integration & Parameter Estimation');   

param_keyword = 'ic_state_vector_apriori';
[ic_state_vector_apriori] = read_param_cfg(orbit_config_fname,param_keyword);
test_ic_apriori = strcmp(ic_state_vector_apriori,'ic');
if test_ic_apriori == 0    
% if test_pod_mode == 1
% IC initialisation based on Short arc solution 
    % N_short_arcs = 0
    N_short_arcs = 3;
    edge_arc_offset = 0;
    % [orbit_arc_long, orbit_arc_short, N_short_arcs] = short_arc_meth (orbit_config_fname, N_short_arcs, -1, 0, 0, edge_arc_offset);
    [orbit_arc_long, orbit_arc_short, N_short_arcs, orbit_model_struct] = short_arc_meth_v2 (orbit_config_fname, N_short_arcs, -1, 0, 0, edge_arc_offset, orbit_model_struct);
    % Overall iterations number
    Niter_estim = Niter_estim + N_short_arcs;    
elseif test_ic_apriori == 1    
    N_short_arcs = 0;
    edge_arc_offset = 0;  
    % [orbit_arc_long, orbit_arc_short, N_short_arcs] = short_arc_meth (orbit_config_fname, N_short_arcs, -1, 0, 0, edge_arc_offset);
    [orbit_arc_long, orbit_arc_short, N_short_arcs, orbit_model_struct] = short_arc_meth_v2 (orbit_config_fname, N_short_arcs, -1, 0, 0, edge_arc_offset, orbit_model_struct);
end

% Orbit Determination iterations loop: start
for iveq = 0 : Niter_estim + 1 
    %fprintf('\n%s %d \n', 'Orbit Integration solution No.: ',iveq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orbit arc length condition
    % [orbit_arc_long, orbit_arc_short, N_short_arcs] = short_arc_meth (orbit_config_fname, N_short_arcs, iveq, orbit_arc_long, orbit_arc_short, edge_arc_offset);
    [orbit_arc_long, orbit_arc_short, N_short_arcs, orbit_model_struct] = short_arc_meth_v2 (orbit_config_fname, N_short_arcs, iveq, orbit_arc_long, orbit_arc_short, edge_arc_offset, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variational Equations numerical integration solution
if iveq < Niter_estim + 1       
% if iveq <= N_short_arcs + 1   % VEQ set to fixed matrices approach      
    % to_VEQ = tic;
    MODEid = sprintf('%s%d','VEQ',iveq);
    VEQ_sol = 1;
    [orbc,orbk,orbt,veqZarray,veqParray,forces_accel_veq] = orbit_integr(orbit_config_fname, VEQ_sol, orbit_model_struct);
    % fprintf('%s %.3f \n', 'Time (min):  VEQ integration:',toc(to_VEQ)/60);
end

% Equation of Motion numerical integration solution
    % to_EQM = tic;
    MODEid = sprintf('%s%d','ESM',iveq);
    VEQ_sol = 0;
    [orbc,orbk,orbt,veqZarray_0,veqParray_0,forces_accel] = orbit_integr(orbit_config_fname, VEQ_sol, orbit_model_struct);
    % fprintf('%s %.3f \n', 'Time (min):  EQM integration:',toc(to_EQM)/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

to_obs_estim = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics / Residuals / Write Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param_estim_01 > 0
    if iveq == Niter_estim + 1 
        i_write = iveq;
    end
    i_write = -1; % iveq;   
    % Observations residuals
    MODEid2 = sprintf('%s%s',MODEid,'obs');
    [rms_orbital,rms_orbc_obs,rms_orbk,rms_orbt,dstn,dorbc_obs,dkepl,dorbt,rms_3D] = mainf_statistout2(GM_glob,orbc,orbk,orbt,obsorbc,obsorbt,veqZarray,veqParray,MODEid2,i_write);
    % fprintf('%s %11.6f %11.6f %11.6f', 'Orbit residuals: RMS(XYZ): ',rms_orbc_obs(1:3));
    % fprintf('\n\n');
    % External Observations residuals (OBS part not used within estimator)
    %MODEid2 = sprintf('%s%s',MODEid,'obsext');
    %[rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,rms_3D] = mainf_statistout2(GM_glob,orbc,orbk,orbt,obsorbc_ext,obsorbt_ext,veqZarray,veqParray,MODEid2,i_write);
    %fprintf('%s %11.6f %11.6f %11.6f \n\n', 'Orbit residuals: XYZ: ',rms_orbc(1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Weights matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
sigma_obs_weight = 1;
% if (iveq > 0) 
%     obs_residuals = dorbc_obs;
%     sigma_obs_weight = [obs_residuals(:,2); obs_residuals(:,3); obs_residuals(:,4)];    
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Observations screening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
test = strcmp(obs_outliers_yn,'y');
if test == 1
    %sigma_obs_weight = pseudo_obs_sigma_outlier;
if (iveq > 0) && (iveq == Niter_estim)
    obs_sigma_test = pseudo_obs_sigma_outlier;
    orbc_xyz = orbc(:,1:4);
    [obs_reduced, obs_removed, Nepochs_removed, obs_residuals, rms_obs] = obs_pseudo_scan(obsorbc, orbc_xyz, obs_sigma_test);
    [a1 a2] = size(orbc_xyz);
    % Gradual incease of threshold to meet maximum limit of removed observations
    Nepochs_removed_limit = (5/100) * a1;
    OBS_limit_step = 1 * 10^-2;    
    while Nepochs_removed > Nepochs_removed_limit
        obs_sigma_test = obs_sigma_test + OBS_limit_step;
        [obs_reduced, obs_removed, Nepochs_removed, obs_residuals, rms_obs] = obs_pseudo_scan(obsorbc, orbc_xyz, obs_sigma_test);
    end
    obsorbc = obs_reduced;
    
% Weights matrix
    % [residuals_d1_pod residuals_d2_pod] = size(obs_residuals);
    % sigma_obs_weight = [obs_residuals(:,2); obs_residuals(:,3); obs_residuals(:,4)];    
    sigma_obs_weight = 1;    
    fprintf('%s %d', 'Observations screening: Nepochs removed: ',Nepochs_removed);
    fprintf('\n\n');
    % fprintf('%s %11.6f %11.6f %11.6f', 'Orbit residuals post obs screening: RMS(XYZ): ',rms_obs(1:3));
    % fprintf('\n\n');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iveq < Niter_estim + 1       
    % Least-squares estimation method
    [Xmatrix,Amatrix2,Wmatrix2] = orbit_estim(orbc, veqZarray, veqParray, obsorbc, sigma_obs_weight, COVPform);
% Estimated Parameters: Update parameters' values to aposteriori values
    ic_apriori_01 = 1;
    [Zo_estim, Xaposteriori, orbit_model_struct] = param_aposteriori_apriori(Xmatrix, ic_apriori_01, orbit_model_struct);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fprintf('%s %.3f \n', 'Time (min): Statistics, OBS screening, Parameters estimation:',toc(to_obs_estim)/60);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Final Observations matrix to output argument
OBS_matrix = obsorbc;
OBS_residuals = dorbc_obs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

% Orbit Determination iterations loop: end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% POD mode : END
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Mode: Orbit Propagation of EQM and VEQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
test = strcmp(orbit_mode,'orbit_propagation_veq');
if test == 1
    %fprintf('\n\n%s\n','Mode: Orbit Propagation EQM and VEQ');     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variational Equations numerical integration solution
    % to_VEQ = tic;
    MODEid = sprintf('%s%d','VEQ',0);
    VEQ_sol = 1;
    [orbc,orbk,orbt,veqZarray,veqParray,forces_accel_veq, Gmatrix, Rmatrix] = orbit_integr(orbit_config_fname, VEQ_sol, orbit_model_struct);
    % fprintf('%s %.3f \n', 'Time (min):  VEQ integration:',toc(to_VEQ)/60);
% Equation of Motion numerical integration solution
    % to_EQM = tic;
    MODEid = sprintf('%s%d','ESM',0);
    VEQ_sol = 0;
    [orbc,orbk,orbt,veqZarray_0,veqParray_0,forces_accel] = orbit_integr(orbit_config_fname, VEQ_sol, orbit_model_struct);
    % fprintf('%s %.3f \n', 'Time (min):  EQM integration:',toc(to_EQM)/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Gmatrix = zeros(3,3); 
    Rmatrix = zeros(3,3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Mode: Orbit Propagation of EQM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
test = strcmp(orbit_mode,'orbit_propagation_eqm');
if test == 1
    %fprintf('\n\n%s\n','Mode: Orbit Propagation EQM');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Equation of Motion numerical integration solution
    % to_EQM = tic;
    MODEid = sprintf('%s%d','ESM',0);
    VEQ_sol = 0;
    [orbc,orbk,orbt,veqZarray_0,veqParray_0,forces_accel] = orbit_integr(orbit_config_fname, VEQ_sol, orbit_model_struct);
    % fprintf('%s %.3f \n', 'Time (min):  EQM integration:',toc(to_EQM)/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    veqZarray = veqZarray_0;
    veqParray = veqParray_0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param_estim_01 == 0     
Xmatrix = 0;
rms_orbc_obs = [0 0 0];
OBS_matrix = 0;
Xaposteriori = 0;
OBS_residuals = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,extorb_crf,extorb_trf] = orbit_ext(orbc,orbk,orbt,veqZarray,veqParray,MODEid,orbit_config_fname, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('%s \n', 'External Orbit Comparison:');
% fprintf('%s ', 'Orbit residuals: RMS(RTN): ' );
% fprintf('%11.6f ', rms_orbital);
% fprintf('\n');
% fprintf('%s ', 'Orbit residuals: RMS(XYZ): ' );
% fprintf('%11.6f ', rms_orbc(1:3));
% fprintf('\n\n')
