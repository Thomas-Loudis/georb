function [gravitational_effects_01, non_gravitational_forces_01, tidal_effects_01, Gravity_field_terms] = prm_dynamics(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_orbdynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Configuration file read : Orbital Dynamics list included/excluded 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Orbit configuration file name 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                                7 July 2022
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational Dynamics   
    param_keyword = 'Earth_Gravity_Field';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    gravitational_effects_yn(1,1) = param_value; % sscanf(param_value,'%d');

    param_keyword = '3rd_body_peturbations';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    gravitational_effects_yn(2,1) = param_value; 
    
    param_keyword = 'Tides_effects';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    gravitational_effects_yn(3,1) = param_value; 
    
    param_keyword = 'Relativity_effects';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    gravitational_effects_yn(4,1) = param_value; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravity field terms 
    param_keyword = 'Gravity_field_terms';
    [Gravity_field_terms] = read_param_cfg(cfg_fname,param_keyword);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tidal effects
    param_keyword = 'solid_earth_tides_1_non_freq';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    tidal_effects_yn(1,1) = param_value; 
    
    param_keyword = 'solid_earth_tides_2_freq';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    tidal_effects_yn(2,1) = param_value; 
 
    param_keyword = 'ocean_tides';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    tidal_effects_yn(3,1) = param_value; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativity effects
    param_keyword = 'Schwarzschild_effect';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    relativistic_effects_yn(1,1) = param_value; 

    param_keyword = 'Lense_Thirring_effect';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    relativistic_effects_yn(2,1) = param_value; 
    
    param_keyword = 'geodesic_effect';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    relativistic_effects_yn(3,1) = param_value; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-Gravitational effects    
    param_keyword = 'non_gravitational_forces';
    [non_gravitational_forces_yn] = read_param_cfg(cfg_fname,param_keyword);
    
% empirical_forces         y   
    param_keyword = 'empirical_forces';
    [empirical_forces_yn] = read_param_cfg(cfg_fname,param_keyword);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .in configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational Dynamics    
    test = strcmp(str1,'Earth_Gravity_Field');
    if test == 1
      gravitational_effects_yn(1,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'3rd_body_peturbations');
    if test == 1
      gravitational_effects_yn(2,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Tides_effects');
    if test == 1
      gravitational_effects_yn(3,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Relativity_effects');
    if test == 1
      gravitational_effects_yn(4,1) = sscanf(line_ith,'%*s %s %*');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravity field terms 
    test = strcmp(str1,'Gravity_field_terms');
    if test == 1
      Gravity_field_terms = sscanf(line_ith,'%*s %s %*');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tidal effects
    test = strcmp(str1,'solid_earth_tides_1_non_freq');
    if test == 1
      tidal_effects_yn(1,1) = sscanf(line_ith,'%*s %s %*');
    end
 
    test = strcmp(str1,'solid_earth_tides_2_freq');
    if test == 1
      tidal_effects_yn(2,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'ocean_tides');
    if test == 1
      tidal_effects_yn(3,1) = sscanf(line_ith,'%*s %s %*');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativity effects
    test = strcmp(str1,'Schwarzschild_effect');
    if test == 1
      relativistic_effects_yn(1,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Lense_Thirring_effect');
    if test == 1
      relativistic_effects_yn(2,1) = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'geodesic_effect');
    if test == 1
      relativistic_effects_yn(3,1) = sscanf(line_ith,'%*s %s %*');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-Gravitational effects    
    test = strcmp(str1,'non_gravitational_forces');
    if test == 1
      non_gravitational_forces_yn = sscanf(line_ith,'%*s %s %*');
    end
    
% empirical_forces         y   
    test = strcmp(str1,'empirical_forces');
    if test == 1
      empirical_forces_yn = sscanf(line_ith,'%*s %s %*');
    end
    
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitational effects list
gravitational_effects_01 = zeros(4,1);
if gravitational_effects_yn(1,1) == 'y'
    gravitational_effects_01(1,1) = 1;
end
if gravitational_effects_yn(2,1) == 'y'
    gravitational_effects_01(2,1) = 1;
end
if gravitational_effects_yn(3,1) == 'y'
    gravitational_effects_01(3,1) = 1;
end
if gravitational_effects_yn(4,1) == 'y'
    gravitational_effects_01(4,1) = 1;
end

% Non-gravitational effects
non_gravitational_forces_01 = 0;
if non_gravitational_forces_yn == 'y'
    non_gravitational_forces_01 = 1;
end
    
% Tidal effects
tidal_effects_01 = zeros(3,1);
if tidal_effects_yn(1,1) == 'y'
    tidal_effects_01(1,1) = 1;
end
if tidal_effects_yn(2,1) == 'y'
    tidal_effects_01(2,1) = 1;
end
if tidal_effects_yn(3,1) == 'y'
    tidal_effects_01(3,1) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
