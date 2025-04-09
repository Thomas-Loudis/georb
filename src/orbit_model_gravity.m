function [orbit_model_matrix] = orbit_model_gravity (orbit_config_matrix,orbit_model_matrix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model_gravity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Parameters/Variables for gravity model and gravity field determination 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               7 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Gravity Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity model data file (spherical harmonic coefficients)
[GM,ae,Cnm,Snm,sCnm,sSnm, n_max_eqm, m_max_eqm, n_max_veq, m_max_veq, tide_system, gfm_struct] = prm_grav_model(orbit_config_matrix);
gfm_struct_glob = gfm_struct;
% GM_glob = GM;
N_param_GRAV = gfm_struct.parameters_number;
orbit_model_matrix.gravity_field = gfm_struct;
orbit_model_matrix.GM_Earth = GM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%--------------------------------------------------------------------------
param_keyword = 'aod_effects';
[aod_effects_yn] = read_param_cfg(orbit_config_matrix,param_keyword);
test = strcmp(aod_effects_yn,'y');
if test == 1 
% Read AOD data
[aod_struct] = prm_aod_data(orbit_config_matrix);
aod_nmax = aod_struct.degree;
aod_struct.degree_order = [aod_nmax aod_nmax];
aod_struct.degree_eqm = [aod_nmax aod_nmax];
aod_struct.degree_veq = [n_max_veq m_max_veq];
orbit_model_matrix.aod_effects = aod_struct;
end
%--------------------------------------------------------------------------
