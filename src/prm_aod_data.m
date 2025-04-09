function [aod_struct] = prm_aod_data(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  prm_aod_data.m
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


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2

% AOD effects y/n
    param_keyword = 'aod_effects';
    [aod_effects_yn] = read_param_cfg(cfg_fname,param_keyword);
    
% AOD effect ID
    param_keyword = 'aod_effect_data_type';
    [aod_effect_ID] = read_param_cfg(cfg_fname,param_keyword);

% AOD data file name  
    param_keyword = 'AOD_data_filename';
    [AOD_data_filename] = read_param_cfg(cfg_fname,param_keyword);

% AOD maximum degree/order  
    param_keyword = 'AOD_degree_max';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    AOD_degree_max_config = sscanf(param_value,'%d %*');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
%cfg_fname
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

% AOD effects y/n
    test = strcmp(str1,'aod_effects');
    if test == 1
      aod_effects_yn = sscanf(line_ith,'%*s %s %*'); 
    end
    
% AOD effect ID
    test = strcmp(str1,'aod_effect_data_type');
    if test == 1
      aod_effect_ID = sscanf(line_ith,'%*s %s %*') ;
    end

% AOD data file name  
    test = strcmp(str1,'AOD_data_filename');
    if test == 1
      AOD_data_filename = sscanf(line_ith,'%*s %s %*')
    end

% AOD maximum degree/order  
    test = strcmp(str1,'AOD_degree_max');
    if test == 1
      AOD_degree_max_config = sscanf(line_ith,'%*s %d %*') ;
    end
    
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(aod_effects_yn,'y');

if test == 1
    
% AOD data read and store to Stokes coefficients matrices
[aod_GM,aod_ae,aod_model_nmax, aod_Cnm, aod_Snm, aod_atm_Cnm, aod_atm_Snm, aod_ocn_Cnm, aod_ocn_Snm, aod_glo_Cnm, aod_glo_Snm, aod_oba_Cnm, aod_oba_Snm, aod_struct] = aod_read(AOD_data_filename);

test_data_effect = strcmp(aod_effect_ID,'ATM');
if test_data_effect == 1    
    Cnm_AOD = aod_atm_Cnm;
    Snm_AOD = aod_atm_Snm;
end

test_data_effect = strcmp(aod_effect_ID,'OCN');
if test_data_effect == 1    
    Cnm_AOD = aod_ocn_Cnm;
    Snm_AOD = aod_ocn_Snm;
end

test_data_effect = strcmp(aod_effect_ID,'GLO');
if test_data_effect == 1    
    Cnm_AOD = aod_glo_Cnm;
    Snm_AOD = aod_glo_Snm;
end

test_data_effect = strcmp(aod_effect_ID,'OBA');
if test_data_effect == 1    
    Cnm_AOD = aod_oba_Cnm;
    Snm_AOD = aod_oba_Snm;
end

if AOD_degree_max_config == -1
    aod_nmax = aod_model_nmax;
else
    aod_nmax = AOD_degree_max_config;
end

else
    aod_GM = 0;
    aod_ae = 0;
    aod_nmax = 0;   
    Cnm_AOD = 0;
    Snm_AOD = 0;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strucutre array
aod_struct.status_yn = aod_effects_yn;
aod_struct.effect_id = aod_effect_ID;
aod_struct.GM     = aod_GM;
aod_struct.radius = aod_ae;
aod_struct.degree = aod_nmax;
aod_struct.Cnm    = Cnm_AOD;
aod_struct.Snm    = Snm_AOD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

