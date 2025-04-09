function [shc_struct] = prm_ocean_tides(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_ocean_tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Data reading and preprocessing: Ocean Tides model file and form
%  spherical harmonic coefficients matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                    27 May 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 07/06/2022  Thomas Loudis Papanikolaou
%             Code minor modifications
% 07/07/2022  Thomas Loudis Papanikolaou
%             Code minor modifications
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2
% Ocean Tides
    param_keyword = 'ocean_tides';
    [ocean_tides_yn] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'ocean_tides_model_fname';
    [ocean_tides_model_fname] = read_param_cfg(cfg_fname,param_keyword);

% Equation of Motion: truncation degree and order
    param_keyword = 'ocean_tides_degree';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    n_max_otides = sscanf(param_value,'%d %*');

    param_keyword = 'ocean_tides_order';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    m_max_otides = sscanf(param_value,'%d %*');
    
% Variational Equations: truncation degree and order
    param_keyword = 'veq_ocean_tides_degree';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    n_max_otides_veq = sscanf(param_value,'%d %*');

    param_keyword = 'veq_ocean_tides_order';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    m_max_otides_veq = sscanf(param_value,'%d %*');    
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
    
% Ocean Tides
    test = strcmp(str1,'ocean_tides');
    if test == 1
      ocean_tides_yn = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'ocean_tides_model_fname');
    if test == 1
      ocean_tides_model_fname = sscanf(line_ith,'%*s %s %*');
    end

% Equation of Motion: truncation degree and order
    test = strcmp(str1,'ocean_tides_degree');
    if test == 1
      n_max_otides = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'ocean_tides_order');
    if test == 1
      m_max_otides = sscanf(line_ith,'%*s %d %*');
    end
    
% Variational Equations: truncation degree and order
    test = strcmp(str1,'veq_ocean_tides_degree');
    if test == 1
      n_max_otides_veq = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'veq_ocean_tides_order');
    if test == 1
      m_max_otides_veq = sscanf(line_ith,'%*s %d %*');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models data read and preprocessing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides model in form of geopotential spherical harmonic coefficients
[shc_struct, delaunay_doodson_multipliers,otides_dCnm_plus,otides_dSnm_plus,otides_dCnm_minus,otides_dSnm_minus] = read_oceantides(ocean_tides_model_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
otides_DelaunayNf = delaunay_doodson_multipliers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
shc_struct.status_yn = ocean_tides_yn;
shc_struct.degree_eqm   = [n_max_otides m_max_otides];
shc_struct.degree_veq   = [n_max_otides_veq m_max_otides_veq];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
