function [Relativity_yn, Relativistic_effects_yn, PPN_beta_parameter, PPN_gama_parameter, cspeedlight] = prm_relativity(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_relativity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Relativistic Effects: Read orbit configuration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Orbit confiugration file to structure
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                    27 May 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2

param_keyword = 'Relativity_effects';
[Relativity_yn] = read_param_cfg(cfg_fname,param_keyword);
    
% Relativistic effects
    param_keyword = 'Schwarzschild_effect';
    [Schwarzschild_effect_yn] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'Lense_Thirring_effect';
    [Lense_Thirring_effect_yn] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'geodesic_effect';
    [geodesic_effect_yn] = read_param_cfg(cfg_fname,param_keyword);

    
    param_keyword = 'PPN_beta_parameter';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PPN_beta_parameter = sscanf(param_value,'%d %*');
          
    param_keyword = 'PPN_gama_parameter';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    PPN_gama_parameter = sscanf(param_value,'%d %*');

    param_keyword = 'C_speed_of_light';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    cspeedlight = sscanf(param_value,'%d %*');    
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
    
% Relativity effects
    test = strcmp(str1,'Schwarzschild_effect');
    if test == 1
      Schwarzschild_effect_yn = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'Lense_Thirring_effect');
    if test == 1
      Lense_Thirring_effect_yn = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'geodesic_effect');
    if test == 1
      geodesic_effect_yn = sscanf(line_ith,'%*s %s %*');
    end

          
    test = strcmp(str1,'PPN_beta_parameter');
    if test == 1
      PPN_beta_parameter = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'PPN_gama_parameter');
    if test == 1
      PPN_gama_parameter = sscanf(line_ith,'%*s %d %*');
    end
    
    test = strcmp(str1,'C_speed_of_light');
    if test == 1
      cspeedlight = sscanf(line_ith,'%*s %d %*');
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Relativistic effects array 
Relativistic_effects_yn(1) = Schwarzschild_effect_yn;
Relativistic_effects_yn(2) = Lense_Thirring_effect_yn;
Relativistic_effects_yn(3) = geodesic_effect_yn;

% test = strcmp(Schwarzschild_effect_yn,'y');
% Relativistic_effects_yn(1,1) = test;
% 
% test = strcmp(Lense_Thirring_effect_yn,'y');
% Relativistic_effects_yn(1,2) = test;
% 
% test = strcmp(geodesic_effect_yn,'y');
% Relativistic_effects_yn(1,3) = test;
