function [orbc,orbk,orbt,veqZarray,veqParray,forces_accel, Gmatrix, Rmatrix] = orbit_integr(cfg_fname, VEQ_sol, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_integr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Main function for the orbit integration for solving the Equation of
%  Motion and the Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 15/04/2021, Thomas Papanikolaou 
%             Function renamed & upgaded from mainf_DOD to orbintegr
% 08/07/2022, Thomas Loudis Papanikolaou 
%             Code modifications for the release of July 2022
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC_MJDo = orbit_model_struct.IC_MJD; 
Zo      = orbit_model_struct.IC_CRF; 
arc     = orbit_model_struct.orbit_arc_length_sec; 
eopdat  = orbit_model_struct.EOP_data; 
dpint   = orbit_model_struct.EOP_interp_no; 
Nparam_GLOB = orbit_model_struct.forces_param_estim_no; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read configuration file
[integr_method, integr_stepsize, integr_order, integr_start, integr_RKN_lamda, integr_RKN_sigma] = prm_integrator(cfg_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrator method ID
test = strcmp(integr_method,'RKN-768');
if test == 1
    intg_method_id = 2;
end

test = strcmp(integr_method,'RKN-654');
if test == 1
    intg_method_id = 3;
end

test = strcmp(integr_method,'Adams-Bashforth');
if test == 1
    intg_method_id = 4;
end

test = strcmp(integr_method,'Adams-Bashforth-Moulton');
if test == 1
    intg_method_id = 5;
end

test = strcmp(integr_method,'Gauss-Jackson');
if test == 1
    intg_method_id = 8;
end

test = strcmp(integr_method,'Gauss-Jackson-pece');
if test == 1
    intg_method_id = 9;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Integrator ID
test = strcmp(integr_start,'RKN-768');
if test == 1
    integr_start_id = 2;
end

test = strcmp(integr_start,'RKN-654');
if test == 1
    integr_start_id = 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-step methods matrix
MSprm(1,1) = intg_method_id;
MSprm(2,1) = integr_order;
MSprm(3,1) = integr_stepsize;
MSprm(4,1) = integr_start_id;

% RKN methods matrix
RKprm(1,1) = integr_start_id;
RKprm(2,1) = integr_stepsize;
RKprm(3,1) = integr_RKN_lamda;
RKprm(4,1) = integr_RKN_sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration of Equation of Motion and Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intg = intg_method_id;
[orbc,err,veqZarray,veqParray,orbk,orbt,forces_accel, Gmatrix, Rmatrix] = intg_main(intg,Zo,arc,MSprm,RKprm,eopdat,dpint,VEQ_sol, orbit_model_struct);
if Nparam_GLOB == 0   
    veqParray = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
