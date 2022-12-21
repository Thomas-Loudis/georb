function [integr_method, integr_stepsize, integr_order, integr_start, integr_RKN_lamda, integr_RKN_sigma] = prm_integrator(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Configuration file read : 
%  Initial Conditions (State Vector and Initial epoch) and orbit arc length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Orbit configuration file name 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment: This function replaced the function prm3 (01/2011) according to
% the new configuration file format for the orbit modelling
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
% Numerical Integration methods 
    % Numerical Integrator
    param_keyword = 'Integration_method';
    [integr_method] = read_param_cfg(cfg_fname,param_keyword);

    % Integration step size
    param_keyword = 'Stepsize';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_stepsize = sscanf(param_value,'%d %*');

% Parameters for multistep methods only    
    param_keyword = 'integrator_order_multistep';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_order = sscanf(param_value,'%d %*');
    
    param_keyword = 'Start_integrator';
    [integr_start] = read_param_cfg(cfg_fname,param_keyword);

% Parameters for Runge-Kutta-Nystrom methods only    
    param_keyword = 'RKN_lamda_coefficient';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_RKN_lamda = sscanf(param_value,'%f %*');

    param_keyword = 'RKN_interpolation_sigma';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    integr_RKN_sigma = sscanf(param_value,'%f %*');    
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
    
% Numerical Integration methods 
    % Numerical Integrator
    test = strcmp(str1,'Integration_method');
    if test == 1
      integr_method = sscanf(line_ith,'%*s %s %*');
    end
    
    test = strcmp(str1,'Stepsize');
    if test == 1
      integr_stepsize = sscanf(line_ith,'%*s %d %*');
    end

% Parameters for multistep methods only    
    test = strcmp(str1,'integrator_order_multistep');
    if test == 1
      integr_order = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'Start_integrator');
    if test == 1
      integr_start = sscanf(line_ith,'%*s %s %*');
    end

% Parameters for Runge-Kutta-Nystrom methods only    
    test = strcmp(str1,'RKN_lamda_coefficient');
    if test == 1
      integr_RKN_lamda = sscanf(line_ith,'%*s %f %*');
    end

    test = strcmp(str1,'RKN_interpolation_sigma');
    if test == 1
      integr_RKN_sigma = sscanf(line_ith,'%*s %f %*');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end
