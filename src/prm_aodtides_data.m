function [aod_tides_struct] = prm_aodtides_data(config_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  prm_aodtides_data.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
% Read Atmosphere and Ocean De-Aliasing (AOD) effects data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname : Orbit modelling configuration file file name
%
% Output arguments:
% - Cnm     : Spherical Harmonic Coefficients (SHC) array
% - Snm     : Spherical Harmonic Coefficients (SHC) array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             3 October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AOD effects y/n
    param_keyword = 'atm_tides';
    [atm_tides_yn] = read_param_cfg(config_struct,param_keyword);

% ATM Data level name 
    param_keyword = 'atm_tides_data_level';
    [atm_tides_data_level] = read_param_cfg(config_struct,param_keyword);

% ATM Data Release 
    param_keyword = 'atm_tides_data_release';
    [atm_tides_data_release] = read_param_cfg(config_struct,param_keyword);
    
% AOD maximum degree/order  
    param_keyword = 'atm_tides_degree';
    [param_value] = read_param_cfg(config_struct,param_keyword);
    atm_tides_degree = sscanf(param_value,'%d %*');        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tides_names_array = [
'P1 163.555 ' 
'K1 165.555 '
'S1 164.555 '
'N2 245.655 ' 
'M2 255.555 '
'L2 265.455 '
'T2 272.556 '
'S2 273.555 ' 
'R2 274.554 ' 
'T3 381.555 '
'S3 382.555 '
'R3 383.555 ' 
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test = strcmp(atm_tides_yn,'y');
if test == 1
%--------------------------------------------------------------------------
% Atmospheric Tides
%--------------------------------------------------------------------------
[k1 k2] = size(tides_names_array);
N_freq = k1;

for i_freq = 1 : N_freq
    
tide_darwin_symbol  = sscanf(tides_names_array(i_freq,:),'%s %*');
doodson_number_char = sscanf(tides_names_array(i_freq,:),'%*s %s %*');

% Data File name according to AOD1B naming conventions
%data_level_id = 'AOD1B';
data_level_id = atm_tides_data_level;
tides_category_id = 'ATM';
%data_release_version = '06';
data_release_version = atm_tides_data_release;
data_prefix = 'asc';
tides_filename = sprintf('%s%s%s%s%s%s%s%s%s', data_level_id,'_', tides_category_id,'_', tide_darwin_symbol,'_', data_release_version,'.',data_prefix);

% Read AOD1B Tides data file
[aod_GM,aod_ae,aod_nmax, dCnm_p, dSnm_p, dCnm_m, dSnm_m] = read_aod_tides(tides_filename);

if i_freq == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
Nmax_tides  = aod_nmax;
Nfreq_tides = N_freq;
delaunay_doodson_multipliers = zeros(Nfreq_tides, 12);
aod_tides_dCnm_plus  = zeros(Nmax_tides+1,Nmax_tides+1,Nfreq_tides);
aod_tides_dSnm_plus  = zeros(Nmax_tides+1,Nmax_tides+1,Nfreq_tides);
aod_tides_dCnm_minus = zeros(Nmax_tides+1,Nmax_tides+1,Nfreq_tides);
aod_tides_dSnm_minus = zeros(Nmax_tides+1,Nmax_tides+1,Nfreq_tides);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Delaunay and Doodson multipliers based on Doodson number
[doodson_multipliers, delaunay_multipliers] = doodson_number(doodson_number_char);

% Array of multipliers
doodson_number_num = sscanf(doodson_number_char,'%f%*');
delaunay_doodson_multipliers (i_freq,:) = [doodson_number_num delaunay_multipliers doodson_multipliers];

% AOD Tides matrices
aod_tides_dCnm_plus (:,:,i_freq) = dCnm_p;
aod_tides_dSnm_plus (:,:,i_freq) = dSnm_p;
aod_tides_dCnm_minus(:,:,i_freq) = dCnm_m;
aod_tides_dSnm_minus(:,:,i_freq) = dSnm_m;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric Tides model from data file in form of geopotential spherical harmonic coefficients
prm_aodtides_data_atm_tides = 0;
if prm_aodtides_data_atm_tides == 2
atm_tides_potential_fname = 'atmosTides_AOD1BRL06.potential.iers.txt';
[atmtides_struct, delaunay_doodson_multipliers,aod_tides_dCnm_plus,aod_tides_dSnm_plus,aod_tides_dCnm_minus,aod_tides_dSnm_minus] = read_oceantides(atm_tides_potential_fname);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%--------------------------------------------------------------------------
else
aod_GM = 0;
aod_ae = 0;
atm_tides_degree = 0;
delaunay_doodson_multipliers = 0;
aod_tides_dCnm_plus = 0;
aod_tides_dSnm_plus = 0;
aod_tides_dCnm_minus = 0;
aod_tides_dSnm_minus = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strucutre array
aod_tides_struct.status_yn                = atm_tides_yn;
aod_tides_struct.GM                       = aod_GM;
aod_tides_struct.radius                   = aod_ae;
aod_tides_struct.degree                   = atm_tides_degree;
aod_tides_struct.delaunay_doodson_multipl = delaunay_doodson_multipliers;
aod_tides_struct.dCnm_plus                = aod_tides_dCnm_plus;
aod_tides_struct.dSnm_plus                = aod_tides_dSnm_plus;
aod_tides_struct.dCnm_minus               = aod_tides_dCnm_minus;
aod_tides_struct.dSnm_minus               = aod_tides_dSnm_minus;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
