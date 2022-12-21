function [emp_cpr_struct] = prm_empirical_cpr (cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_empirical_cpr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read configuration file for Empirical Forces based on cycle-per-revolution terms (periodic accelerations) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                    6 July 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
% 14/12/2022  Dr. Thomas Loudis Papanikolaou
%             Code upgrade to call function empirical_cpr_init and pass the output arguments to strucutre arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2

    % Empirical Forces modelling
    param_keyword = 'empirical_forces';
    [empirical_forces_yn] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'empirical_frame';
    [empirical_frame] = read_param_cfg(cfg_fname,param_keyword);
    
% Empirical Forces' parameters matrix to be considered yes/no
    
% Bias empirical accelerations per axis
    param_keyword = 'empirical_bias_axis1';
    [empirical_bias_axis1] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'empirical_bias_axis2';
    [empirical_bias_axis2] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'empirical_bias_axis3';
    [empirical_bias_axis3] = read_param_cfg(cfg_fname,param_keyword);
    
% Cycle-per-revolution terms: Cosine (C) and Sine (S) coefficients per axis  
    param_keyword = 'cpr_C_axis1';
    [cpr_C_axis1] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'cpr_S_axis1';
    [cpr_S_axis1] = read_param_cfg(cfg_fname,param_keyword);

    
    param_keyword = 'cpr_C_axis2';
    [cpr_C_axis2] = read_param_cfg(cfg_fname,param_keyword);
    
    param_keyword = 'cpr_S_axis2';
    [cpr_S_axis2] = read_param_cfg(cfg_fname,param_keyword);

    
    param_keyword = 'cpr_C_axis3';
    [cpr_C_axis3] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'cpr_S_axis3';
    [cpr_S_axis3] = read_param_cfg(cfg_fname,param_keyword);


% Number of cycles per revolution (once, twice, n per revolution)
    param_keyword = 'cpr_freq_number';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    cpr_freq_number = sscanf(param_value,'%d %*');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

% Empirical Forces modelling
    test = strcmp(str1,'empirical_forces');
    if test == 1
      empirical_forces_yn = sscanf(line_ith,'%*s %s %*') ;
    end
    
    test = strcmp(str1,'empirical_frame');
    if test == 1
      empirical_frame = sscanf(line_ith,'%*s %s %*') ;
    end

% Empirical Forces' parameters matrix to be considered yes/no
    
% Bias empirical accelerations per axis
    test = strcmp(str1,'empirical_bias_axis1');
    if test == 1
      empirical_bias_axis1 = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'empirical_bias_axis2');
    if test == 1
      empirical_bias_axis2 = sscanf(line_ith,'%*s %s %*') ;
    end
    
    test = strcmp(str1,'empirical_bias_axis3');
    if test == 1
      empirical_bias_axis3 = sscanf(line_ith,'%*s %s %*') ;
    end

% Cycle-per-revolution terms: Cosine (C) and Sine (S) coefficients per axis  
    test = strcmp(str1,'cpr_C_axis1');
    if test == 1
      cpr_C_axis1 = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'cpr_S_axis1');
    if test == 1
      cpr_S_axis1 = sscanf(line_ith,'%*s %s %*') ;
    end
    
    
    test = strcmp(str1,'cpr_C_axis2');
    if test == 1
      cpr_C_axis2 = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'cpr_S_axis2');
    if test == 1
      cpr_S_axis2 = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'cpr_C_axis3');
    if test == 1
      cpr_C_axis3 = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'cpr_S_axis3');
    if test == 1
      cpr_S_axis3 = sscanf(line_ith,'%*s %s %*') ;
    end

% Number of cycles per revolution (once, twice, n per revolution)
    test = strcmp(str1,'cpr_freq_number');
    if test == 1
      cpr_freq_number = sscanf(line_ith,'%*s %d %*') ;
    end
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMP_FORCE_PARAM_yn : Array for defining the parameters to be estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format: EMP_FORCE_PARAM_yn(ei, :) = [bias C_nCPR S_nCPR ... ]
%         ei : axis i (1,2,3) e.g. radial, along-track, cross-track
% e1 [bias C_nCPR S_nCPR ... ]
% e2 [bias C_nCPR S_nCPR ... ]
% e3 [bias C_nCPR S_nCPR ... ]

% Preallocation
EMP_FORCE_PARAM_yn   = zeros(3,3);

if empirical_forces_yn == 'y'
    
if empirical_bias_axis1 == 'y'
    EMP_FORCE_PARAM_yn(1,1) = 1;
end
if empirical_bias_axis2 == 'y'
    EMP_FORCE_PARAM_yn(2,1) = 1;
end
if empirical_bias_axis3 == 'y'
    EMP_FORCE_PARAM_yn(3,1) = 1;
end

if cpr_C_axis1 == 'y'
    EMP_FORCE_PARAM_yn(1,2) = 1;
end

if cpr_S_axis1 == 'y'
    EMP_FORCE_PARAM_yn(1,3) = 1;
end


if cpr_C_axis2 == 'y'
    EMP_FORCE_PARAM_yn(2,2) = 1;
end

if cpr_S_axis2 == 'y'
    EMP_FORCE_PARAM_yn(2,3) = 1;
end


if cpr_C_axis3 == 'y'
    EMP_FORCE_PARAM_yn(3,2) = 1;
end

if cpr_S_axis3 == 'y'
    EMP_FORCE_PARAM_yn(3,3) = 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces : y/n effects to be considered
EMP_FORCE_effect_01 = strcmp(empirical_forces_yn,'y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces parameters: Number of unknown parameters to be estimated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d1, d2] = size(EMP_FORCE_PARAM_yn);
Nparam_EMP_FORCE = 0;
for i1 = 1 : d1
    for i2 = 1 : d2
        if EMP_FORCE_PARAM_yn(i1,i2) == 1
            Nparam_EMP_FORCE = Nparam_EMP_FORCE + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of empirical forces model
[accel_empirical] = empirical_cpr_init (EMP_FORCE_PARAM_yn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces CPR: y/n 
emp_cpr_struct.effect_01 = EMP_FORCE_effect_01;
% Array for defining the parameters to be estimated
emp_cpr_struct.parameters_01 = EMP_FORCE_PARAM_yn;
% Reference Frame of the empirical acclerations
emp_cpr_struct.reference_frame = empirical_frame;
% Empirical accelerations vector matrix 3x3 for radial, along-track, cross-track
emp_cpr_struct.acceleration_matrix = accel_empirical;
% Number of parameters of empirical forces
emp_cpr_struct.parameters_number = Nparam_EMP_FORCE;
% Fundamental frequency of the cycle-per-revolution terms
emp_cpr_struct.cpr_frequency_number = cpr_freq_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
